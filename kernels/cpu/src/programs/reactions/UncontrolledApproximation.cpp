/**
 * << detailed description >>
 *
 * @file UncontrolledApproximation.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#include <readdy/kernel/cpu/programs/reactions/UncontrolledApproximation.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

using particle_type = readdy::model::Particle;

UncontrolledApproximation::UncontrolledApproximation(const CPUKernel *const kernel)
        : kernel(kernel) {

}

void UncontrolledApproximation::execute() {
    const auto &ctx = kernel->getKernelContext();
    const auto &fixPos = ctx.getFixPositionFun();
    const auto &dt = ctx.getTimeStep();
    auto data = kernel->getKernelStateModel().getParticleData();
    std::vector<particle_type> newParticles{};
    std::vector<std::function<void()>> events{};

    // reactions with one educt
    {
        auto it_type = data->begin_types();

        while (it_type != data->end_types()) {
            // gather reactions
            const auto &reactions = ctx.getOrder1Reactions(*it_type);
            for (const auto &reaction : reactions) {
                auto r = reaction->getRate() * dt;
                if (readdy::model::rnd::uniform() < r) {
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
                                particle_type outParticle1{};
                                reaction->perform(particle, particle, outParticle1, outParticle1);
                                newParticles.push_back(outParticle1);
                                break;
                            }
                            case 2: {
                                const auto particle = (*_data)[particleIdx];
                                particle_type outParticle1{}, outParticle2{};
                                reaction->perform(particle, particle, outParticle1, outParticle2);
                                newParticles.push_back(outParticle1);
                                newParticles.push_back(outParticle2);
                                break;
                            }
                            default: {
                                log::console()->error("This should not happen!");
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
            for (const auto &neighbor : it->second) {
                if (idx1 > neighbor.idx) continue;
                const auto &reactions = ctx.getOrder2Reactions(
                        *(data->begin_types() + idx1), *(data->begin_types() + neighbor.idx)
                );

                const auto distSquared = neighbor.d2;
                for (const auto &reaction : reactions) {
                    // if close enough and coin flip successful
                    if (distSquared < reaction->getEductDistanceSquared()
                        && readdy::model::rnd::uniform() < reaction->getRate() * dt) {
                        events.push_back([idx1, neighbor, this, &newParticles, &reaction] {
                            auto &&_data = kernel->getKernelStateModel().getParticleData();
                            if (_data->isMarkedForDeactivation(idx1)) return;
                            if (_data->isMarkedForDeactivation(neighbor.idx)) return;
                            _data->markForDeactivation(idx1);
                            _data->markForDeactivation(neighbor.idx);
                            const auto inParticle1 = (*_data)[idx1];
                            const auto inParticle2 = (*_data)[neighbor.idx];
                            switch (reaction->getNProducts()) {
                                case 1: {
                                    particle_type out{};
                                    reaction->perform(inParticle1, inParticle2, out, out);
                                    newParticles.push_back(out);
                                    break;
                                }
                                case 2: {
                                    particle_type out1{}, out2{};
                                    reaction->perform(inParticle1, inParticle2, out1, out2);
                                    newParticles.push_back(out1);
                                    newParticles.push_back(out2);
                                    break;
                                }
                                default: {
                                    log::console()->error("This should not happen!");
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
                  [&fixPos](particle_type &p) { fixPos(p.getPos()); });

    // update data structure
    data->deactivateMarked();
    data->addParticles(newParticles);
}
}
}
}
}
}