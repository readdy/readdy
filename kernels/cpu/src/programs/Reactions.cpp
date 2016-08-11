/**
 * << detailed description >>
 *
 * @file CPUDefaultReactionProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/programs/Reactions.h>
#include <readdy/model/RandomProvider.h>

using _rdy_particle_t = readdy::model::Particle;

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                namespace reactions {
                    UncontrolledApproximation::UncontrolledApproximation(const CPUKernel *const kernel)
                            : kernel(kernel) {

                    }

                    void UncontrolledApproximation::execute() {
                        const auto &ctx = kernel->getKernelContext();
                        const auto &dist = ctx.getDistSquaredFun();
                        const auto &fixPos = ctx.getFixPositionFun();
                        const auto &dt = ctx.getTimeStep();
                        auto data = kernel->getKernelStateModel().getParticleData();
                        auto rnd = readdy::model::RandomProvider();
                        std::vector<_rdy_particle_t> newParticles{};
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
                                                    const auto particle = (*_data)[particleIdx];
                                                    _rdy_particle_t outParticle1{};
                                                    reaction->perform(particle, particle, outParticle1, outParticle1);
                                                    newParticles.push_back(outParticle1);
                                                    break;
                                                }
                                                case 2: {
                                                    const auto particle = (*_data)[particleIdx];
                                                    _rdy_particle_t outParticle1{}, outParticle2{};
                                                    reaction->perform(particle, particle, outParticle1, outParticle2);
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
                            auto cbegin = kernel->getKernelStateModel().getNeighborList()->pairs->cbegin();
                            const auto cend = kernel->getKernelStateModel().getNeighborList()->pairs->cend();
                            for (auto &it = cbegin; it != cend; ++it) {
                                const auto idx1 = it->first;
                                for(const auto& idx2 : it->second) {
                                    if(idx1 > idx2) continue;
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
                                                        _rdy_particle_t out{};
                                                        reaction->perform(inParticle1, inParticle2, out, out);
                                                        newParticles.push_back(out);
                                                        break;
                                                    }
                                                    case 2: {
                                                        _rdy_particle_t out1{}, out2{};
                                                        reaction->perform(inParticle1, inParticle2, out1, out2);
                                                        newParticles.push_back(out1);
                                                        newParticles.push_back(out2);
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
                        }
                        // shuffle reactions
                        std::random_shuffle(events.begin(), events.end());

                        // execute reactions
                        std::for_each(events.begin(), events.end(), [](const std::function<void()> &f) { f(); });

                        // reposition particles to respect the periodic b.c.
                        std::for_each(newParticles.begin(), newParticles.end(),
                                      [&fixPos](_rdy_particle_t &p) { fixPos(p.getPos()); });

                        // update data structure
                        data->deactivateMarked();
                        data->addParticles(newParticles);
                    }
                }

            }
        }
    }
}