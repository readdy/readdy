/**
 * << detailed description >>
 *
 * @file SingleCPUDefaultReactionProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#include <readdy/kernel/singlecpu/programs/SingleCPUDefaultReactionProgram.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {

                struct ReactionEvent {
                    unsigned int order;
                    // TODO: do we need these?
                    unsigned long idx1, idx2;
                    std::function<void()> action;
                };

                SingleCPUDefaultReactionProgram::SingleCPUDefaultReactionProgram(SingleCPUKernel const *const kernel) : kernel(kernel) {

                }

                void SingleCPUDefaultReactionProgram::execute() {
                    const auto &ctx = kernel->getKernelContext();
                    const auto &dt = ctx.getTimeStep();
                    auto data = kernel->getKernelStateModelSingleCPU().getParticleData();
                    auto &rnd = kernel->getRandomProvider();
                    std::vector<readdy::model::Particle> particlesToBeAdded{};
                    std::vector<ReactionEvent> reactionEvents{};

                    // reactions with one educt
                    {
                        auto it_type = data->begin_types();

                        while (it_type != data->end_types()) {
                            // gather reactions
                            const auto &reactions = ctx.getOrder1Reactions(*it_type);

                            for (const auto &reaction : reactions) {
                                if (rnd.getUniform() < reaction->getRate() * dt) {
                                    const unsigned long particleIdx = (const unsigned long) (it_type - data->begin_types());

                                    ReactionEvent evt{};
                                    evt.order = 1;
                                    evt.idx1 = particleIdx;
                                    evt.action = [&] {
                                        if (data->isMarkedForDeactivation(particleIdx)) return;

                                        data->markForDeactivation(particleIdx);
                                        switch (reaction->getNProducts()) {
                                            case 0: {
                                                // no op
                                                break;
                                            }
                                            case 1: {
                                                if (mapping_11.find(reaction->getId()) != mapping_11.end()) {
                                                    particlesToBeAdded.push_back(mapping_11[reaction->getId()]((*data)[particleIdx]));
                                                } else {
                                                    const auto particle = (*data)[particleIdx];
                                                    readdy::model::Particle p1{};
                                                    reaction->perform(particle, particle, p1, p1);
                                                    particlesToBeAdded.push_back(p1);
                                                }
                                                break;
                                            }
                                            case 2: {
                                                const auto particle = (*data)[particleIdx];
                                                readdy::model::Particle p1{}, p2{};
                                                if (mapping_12.find(reaction->getId()) != mapping_12.end()) {
                                                    mapping_12[reaction->getId()](particle, p1, p2);
                                                } else {
                                                    reaction->perform(particle, particle, p1, p2);
                                                }
                                                particlesToBeAdded.push_back(p1);
                                                particlesToBeAdded.push_back(p2);
                                                break;
                                            }
                                            default: {
                                                BOOST_LOG_TRIVIAL(error) << "This should not happen!";
                                            }
                                        }
                                    };
                                    reactionEvents.push_back(std::move(evt));
                                }
                            }
                            ++it_type;
                        }
                    }

                    // reactions with two educts
                    {
                        const auto *neighborList = kernel->getKernelStateModelSingleCPU().getNeighborList();
                        for (auto &&it = neighborList->begin(); it != neighborList->end(); ++it) {
                            const auto idx1 = it->idx1, idx2 = it->idx2;
                            const auto &reactions = ctx.getOrder2Reactions(*(data->begin_types() + idx1), *(data->begin_types() + idx2));
                            // todo: make this generic (periodic bc)
                            const auto posDiff = *(data->begin_positions() + idx1) - *(data->begin_positions() + idx2);
                            const auto distSquared = posDiff * posDiff;

                            for (const auto &reaction : reactions) {
                                // if close enough and coin flip successful
                                if (distSquared < reaction->getEductDistance() * reaction->getEductDistance() && rnd.getUniform() < reaction->getRate() * dt) {
                                    ReactionEvent evt{};
                                    evt.order = 2;
                                    evt.idx1 = idx1;
                                    evt.idx2 = idx2;
                                    evt.action = [&] {
                                        if (data->isMarkedForDeactivation(idx1)) return;
                                        if (data->isMarkedForDeactivation(idx2)) return;
                                        data->markForDeactivation(idx1);
                                        data->markForDeactivation(idx2);
                                        const auto in1 = (*data)[idx1];
                                        const auto in2 = (*data)[idx2];
                                        switch (reaction->getNProducts()) {
                                            case 1: {
                                                if (mapping_21.find(reaction->getId()) != mapping_21.end()) {
                                                    particlesToBeAdded.push_back(mapping_21[reaction->getId()](in1, in2));
                                                } else {
                                                    readdy::model::Particle p1{};
                                                    reaction->perform(in1, in2, p1, p1);
                                                    particlesToBeAdded.push_back(p1);
                                                }
                                                break;
                                            }
                                            case 2: {
                                                readdy::model::Particle p1{}, p2{};
                                                if (mapping_22.find(reaction->getId()) != mapping_22.end()) {
                                                    mapping_22[reaction->getId()](in1, in2, p1, p2);
                                                } else {
                                                    reaction->perform(in1, in2, p1, p2);
                                                }
                                                particlesToBeAdded.push_back(p1);
                                                particlesToBeAdded.push_back(p2);
                                                break;
                                            }
                                            default: {
                                                BOOST_LOG_TRIVIAL(error) << "This should not happen!";
                                            }
                                        }
                                    };

                                    reactionEvents.push_back(std::move(evt));
                                }
                            }
                        }
                    }
                    // shuffle reactions
                    std::random_shuffle(reactionEvents.begin(), reactionEvents.end());

                    // execute reactions
                    for (auto &&event : reactionEvents) {
                        event.action();
                    }
                    // update data structure
                    data->deactivateMarked();
                    data->addParticles(particlesToBeAdded);
                }

                void SingleCPUDefaultReactionProgram::configure() {

                }

                void SingleCPUDefaultReactionProgram::registerReactionScheme_11(const std::string &reactionName, reaction_11 fun) {
                    auto reaction = kernel->getKernelContext().getReactionOrder1WithName(reactionName);
                    if (reaction) {
                        if (reaction->getNEducts() == 1 && reaction->getNProducts() == 1) {
                            mapping_11.emplace(reaction->getId(), fun);
                        } else {
                            throw std::runtime_error("Reaction did not have exactly one product and exactly one educt.");
                        }
                    } else {
                        throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
                    }
                }

                void SingleCPUDefaultReactionProgram::registerReactionScheme_12(const std::string &reactionName, reaction_12 fun) {
                    auto reaction = kernel->getKernelContext().getReactionOrder1WithName(reactionName);
                    if (reaction) {
                        if (reaction->getNEducts() == 1 && reaction->getNProducts() == 2) {
                            mapping_12.emplace(reaction->getId(), fun);
                        } else {
                            throw std::runtime_error("Reaction did not have exactly two products and exactly one educt.");
                        }
                    } else {
                        throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
                    }
                }

                void SingleCPUDefaultReactionProgram::registerReactionScheme_21(const std::string &reactionName, reaction_21 fun) {
                    auto reaction = kernel->getKernelContext().getReactionOrder2WithName(reactionName);
                    if (reaction) {
                        if (reaction->getNEducts() == 2 && reaction->getNProducts() == 1) {
                            mapping_21.emplace(reaction->getId(), fun);
                        } else {
                            throw std::runtime_error("Reaction did not have exactly one product and exactly two educts.");
                        }
                    } else {
                        throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
                    }
                }

                void SingleCPUDefaultReactionProgram::registerReactionScheme_22(const std::string &reactionName, reaction_22 fun) {
                    auto reaction = kernel->getKernelContext().getReactionOrder2WithName(reactionName);
                    if (reaction) {
                        if (reaction->getNEducts() == 2 && reaction->getNProducts() == 2) {
                            mapping_22.emplace(reaction->getId(), fun);
                        } else {
                            throw std::runtime_error("Reaction did not have exactly two products and exactly two educts.");
                        }
                    } else {
                        throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
                    }
                }


            }
        }
    }
}


