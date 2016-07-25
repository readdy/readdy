/**
 * << detailed description >>
 *
 * @file SingleCPUDefaultReactionProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>

using particle_t = readdy::model::Particle;

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {

                namespace reactions {
                    UncontrolledApproximation::UncontrolledApproximation(SingleCPUKernel const *const kernel)
                            : kernel(kernel) {
                    }

                    void UncontrolledApproximation::execute() {
                        const auto &ctx = kernel->getKernelContext();
                        const auto &dist = ctx.getDistSquaredFun();
                        const auto &fixPos = ctx.getFixPositionFun();
                        const auto &dt = ctx.getTimeStep();
                        auto data = kernel->getKernelStateModel().getParticleData();
                        auto &rnd = kernel->getRandomProvider();
                        std::vector<particle_t> newParticles{};
                        std::vector<std::function<void()>> events{};

                        // reactions with one educt
                        {
                            auto it_type = data->begin_types();

                            while (it_type != data->end_types()) {
                                // gather reactions
                                const auto &reactions = ctx.getOrder1Reactions(*it_type);
                                for (const auto &reaction : reactions) {
                                    auto r = reaction->getRate() * dt;
                                    if (rnd.getUniform() < r) {
                                        const size_t particleIdx = (const size_t) (it_type - data->begin_types());
                                        events.push_back([particleIdx, &newParticles, &reaction, this] {
                                            auto &&_data = kernel->getKernelStateModel().getParticleData();
                                            if (_data->isMarkedForDeactivation(particleIdx)) return;
                                            _data->markForDeactivation(particleIdx);

                                            switch (reaction->getNProducts()) {
                                                case 0: {
                                                    // no operation, just deactivation
                                                    // (read out loud with saxony accent)
                                                    break;
                                                }
                                                case 1: {
                                                    if (mapping_11.find(reaction->getId()) != mapping_11.end()) {
                                                        newParticles.push_back(
                                                                mapping_11[reaction->getId()]((*_data)[particleIdx])
                                                        );
                                                    } else {
                                                        const auto particle = (*_data)[particleIdx];
                                                        particle_t outParticle1{};
                                                        reaction->perform(particle, particle, outParticle1,
                                                                          outParticle1);
                                                        newParticles.push_back(outParticle1);
                                                    }
                                                    break;
                                                }
                                                case 2: {
                                                    const auto particle = (*_data)[particleIdx];
                                                    particle_t outParticle1{}, outParticle2{};
                                                    if (mapping_12.find(reaction->getId()) != mapping_12.end()) {
                                                        mapping_12[reaction->getId()](particle, outParticle1,
                                                                                      outParticle2);
                                                    } else {
                                                        reaction->perform(particle, particle, outParticle1,
                                                                          outParticle2);
                                                    }
                                                    newParticles.push_back(outParticle1);
                                                    newParticles.push_back(outParticle2);
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
                            const auto *neighborList = kernel->getKernelStateModel().getNeighborList();
                            for (auto &&it = neighborList->begin(); it != neighborList->end(); ++it) {
                                const auto idx1 = it->idx1, idx2 = it->idx2;
                                const auto &reactions = ctx.getOrder2Reactions(
                                        *(data->begin_types() + idx1), *(data->begin_types() + idx2)
                                );

                                const auto distSquared = dist(
                                        *(data->begin_positions() + idx1), *(data->begin_positions() + idx2)
                                );

                                for (const auto &reaction : reactions) {
                                    // if close enough and coin flip successful
                                    if (distSquared < reaction->getEductDistance() * reaction->getEductDistance()
                                        && rnd.getUniform() < reaction->getRate() * dt) {
                                        events.push_back([idx1, idx2, this, &newParticles, &reaction] {
                                            auto &&_data = kernel->getKernelStateModel().getParticleData();
                                            if (_data->isMarkedForDeactivation(idx1)) return;
                                            if (_data->isMarkedForDeactivation(idx2)) return;
                                            _data->markForDeactivation(idx1);
                                            _data->markForDeactivation(idx2);
                                            const auto inParticle1 = (*_data)[idx1];
                                            const auto inParticle2 = (*_data)[idx2];
                                            switch (reaction->getNProducts()) {
                                                case 1: {
                                                    if (mapping_21.find(reaction->getId()) != mapping_21.end()) {
                                                        newParticles.push_back(
                                                                mapping_21[reaction->getId()](inParticle1,
                                                                                              inParticle2));
                                                    } else {
                                                        particle_t outParticle1{};
                                                        reaction->perform(inParticle1, inParticle2, outParticle1,
                                                                          outParticle1);
                                                        newParticles.push_back(outParticle1);
                                                    }
                                                    break;
                                                }
                                                case 2: {
                                                    particle_t outParticle1{}, outParticle2{};
                                                    if (mapping_22.find(reaction->getId()) != mapping_22.end()) {
                                                        mapping_22[reaction->getId()](inParticle1, inParticle2,
                                                                                      outParticle1, outParticle2);
                                                    } else {
                                                        reaction->perform(inParticle1, inParticle2, outParticle1,
                                                                          outParticle2);
                                                    }
                                                    newParticles.push_back(outParticle1);
                                                    newParticles.push_back(outParticle2);
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
                        std::random_shuffle(events.begin(), events.end());

                        // execute reactions
                        std::for_each(events.begin(), events.end(), [](const std::function<void()> &f) { f(); });

                        // reposition particles to respect the periodic b.c.
                        std::for_each(newParticles.begin(), newParticles.end(),
                                      [&fixPos](particle_t &p) { fixPos(p.getPos()); });

                        // update data structure
                        data->deactivateMarked();
                        data->addParticles(newParticles);
                    }

                    void UncontrolledApproximation::registerReactionScheme_11(const std::string &reactionName,
                                                                              reaction_11 fun) {
                        auto reaction = kernel->getKernelContext().getReactionOrder1WithName(reactionName);
                        if (reaction) {
                            if (reaction->getNEducts() == 1 && reaction->getNProducts() == 1) {
                                mapping_11.emplace(reaction->getId(), fun);
                            } else {
                                throw std::runtime_error(
                                        "Reaction did not have exactly one product and exactly one educt.");
                            }
                        } else {
                            throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
                        }
                    }

                    void UncontrolledApproximation::registerReactionScheme_12(const std::string &reactionName,
                                                                              reaction_12 fun) {
                        auto reaction = kernel->getKernelContext().getReactionOrder1WithName(reactionName);
                        if (reaction) {
                            if (reaction->getNEducts() == 1 && reaction->getNProducts() == 2) {
                                mapping_12.emplace(reaction->getId(), fun);
                            } else {
                                throw std::runtime_error(
                                        "Reaction did not have exactly two products and exactly one educt.");
                            }
                        } else {
                            throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
                        }
                    }

                    void UncontrolledApproximation::registerReactionScheme_21(const std::string &reactionName,
                                                                              reaction_21 fun) {
                        auto reaction = kernel->getKernelContext().getReactionOrder2WithName(reactionName);
                        if (reaction) {
                            if (reaction->getNEducts() == 2 && reaction->getNProducts() == 1) {
                                mapping_21.emplace(reaction->getId(), fun);
                            } else {
                                throw std::runtime_error(
                                        "Reaction did not have exactly one product and exactly two educts.");
                            }
                        } else {
                            throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
                        }
                    }

                    void UncontrolledApproximation::registerReactionScheme_22(const std::string &reactionName,
                                                                              reaction_22 fun) {
                        auto reaction = kernel->getKernelContext().getReactionOrder2WithName(reactionName);
                        if (reaction) {
                            if (reaction->getNEducts() == 2 && reaction->getNProducts() == 2) {
                                mapping_22.emplace(reaction->getId(), fun);
                            } else {
                                throw std::runtime_error(
                                        "Reaction did not have exactly two products and exactly two educts.");
                            }
                        } else {
                            throw std::runtime_error("No reaction with name \"" + reactionName + "\" registered.");
                        }
                    }


                }
            }
        }
    }
}


