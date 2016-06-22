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

                SingleCPUDefaultReactionProgram::SingleCPUDefaultReactionProgram(SingleCPUKernel const *const kernel) : kernel(kernel) {

                }

                void SingleCPUDefaultReactionProgram::execute() {
                    const auto &ctx = kernel->getKernelContext();
                    const auto &dist = ctx.getDistSquaredFun();
                    const auto &fixPos = ctx.getFixPositionFun();
                    const auto &dt = ctx.getTimeStep();
                    auto data = kernel->getKernelStateModelSingleCPU().getParticleData();
                    auto &rnd = kernel->getRandomProvider();
                    std::vector<readdy::model::Particle> particlesToBeAdded{};
                    std::vector<std::function<void()>> reactionEvents{};

                    // reactions with one educt
                    {
                        auto it_type = data->begin_types();

                        while (it_type != data->end_types()) {
                            // gather reactions
                            const auto &reactions = ctx.getOrder1Reactions(*it_type);

                            for (const auto &reaction : reactions) {
                                if (rnd.getUniform() < reaction->getRate() * dt) {
                                    const unsigned long particleIdx = (const unsigned long) (it_type - data->begin_types());

                                    reactionEvents.push_back([&] {
                                        if (data->isMarkedForDeactivation(particleIdx)) return;

                                        data->markForDeactivation(particleIdx);
                                        switch (reaction->getNProducts()) {
                                            case 0: {
                                                // no operation, just deactivation
                                                // (read out loud with saxony accent)
                                                break;
                                            }
                                            case 1: {
                                                if (mapping_11.find(reaction->getId()) != mapping_11.end()) {
                                                    particlesToBeAdded.push_back(mapping_11[reaction->getId()]((*data)[particleIdx]));
                                                } else {
                                                    const auto particle = (*data)[particleIdx];
                                                    readdy::model::Particle outParticle1{};
                                                    reaction->perform(particle, particle, outParticle1, outParticle1);
                                                    particlesToBeAdded.push_back(outParticle1);
                                                }
                                                break;
                                            }
                                            case 2: {
                                                const auto particle = (*data)[particleIdx];
                                                readdy::model::Particle outParticle1{}, outParticle2{};
                                                if (mapping_12.find(reaction->getId()) != mapping_12.end()) {
                                                    mapping_12[reaction->getId()](particle, outParticle1, outParticle2);
                                                } else {
                                                    reaction->perform(particle, particle, outParticle1, outParticle2);
                                                }
                                                particlesToBeAdded.push_back(outParticle1);
                                                particlesToBeAdded.push_back(outParticle2);
                                                break;
                                            }
                                            default: {
                                                BOOST_LOG_TRIVIAL(error) << "This should not happen!";
                                            }
                                        }
                                    });
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

                            const auto distSquared = dist(*(data->begin_positions() + idx1), *(data->begin_positions() + idx2));

                            for (const auto &reaction : reactions) {
                                // if close enough and coin flip successful
                                if (distSquared < reaction->getEductDistance() * reaction->getEductDistance() && rnd.getUniform() < reaction->getRate() * dt) {
                                    reactionEvents.push_back([&] {
                                        if (data->isMarkedForDeactivation(idx1)) return;
                                        if (data->isMarkedForDeactivation(idx2)) return;
                                        data->markForDeactivation(idx1);
                                        data->markForDeactivation(idx2);
                                        const auto inParticle1 = (*data)[idx1];
                                        const auto inParticle2 = (*data)[idx2];
                                        switch (reaction->getNProducts()) {
                                            case 1: {
                                                if (mapping_21.find(reaction->getId()) != mapping_21.end()) {
                                                    particlesToBeAdded.push_back(mapping_21[reaction->getId()](inParticle1, inParticle2));
                                                } else {
                                                    readdy::model::Particle outParticle1{};
                                                    reaction->perform(inParticle1, inParticle2, outParticle1, outParticle1);
                                                    particlesToBeAdded.push_back(outParticle1);
                                                }
                                                break;
                                            }
                                            case 2: {
                                                readdy::model::Particle outParticle1{}, outParticle2{};
                                                if (mapping_22.find(reaction->getId()) != mapping_22.end()) {
                                                    mapping_22[reaction->getId()](inParticle1, inParticle2, outParticle1, outParticle2);
                                                } else {
                                                    reaction->perform(inParticle1, inParticle2, outParticle1, outParticle2);
                                                }
                                                particlesToBeAdded.push_back(outParticle1);
                                                particlesToBeAdded.push_back(outParticle2);
                                                break;
                                            }
                                            default: {
                                                BOOST_LOG_TRIVIAL(error) << "This should not happen!";
                                            }
                                        }
                                    });
                                }
                            }
                        }
                    }
                    // shuffle reactions
                    std::random_shuffle(reactionEvents.begin(), reactionEvents.end());

                    // execute reactions
                    for (auto &&event : reactionEvents) {
                        event();
                    }
                    // reposition particles to respect the periodic b.c.
                    std::for_each(particlesToBeAdded.begin(), particlesToBeAdded.end(), [&fixPos](readdy::model::Particle &p) { fixPos(p.getPos()); });
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


