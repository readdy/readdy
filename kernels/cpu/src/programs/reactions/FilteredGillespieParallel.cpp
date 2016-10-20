#include <readdy/kernel/cpu/programs/reactions/FilteredGillespieParallel.h>
#include <future>

/**
 * << detailed description >>
 *
 * @file FilteredGillespieParallel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

readdy::kernel::cpu::programs::reactions::FilteredGillespieParallel::FilteredGillespieParallel(
        const readdy::kernel::cpu::CPUKernel *const kernel) : GillespieParallel(kernel) {}

void readdy::kernel::cpu::programs::reactions::FilteredGillespieParallel::handleBoxReactions() {
    using promise_t = std::promise<std::set<unsigned long>>;
    using promise_new_particles_t = std::promise<std::vector<particle_t>>;

    auto worker = [this](SlicedBox &box, ctx_t ctx, data_t data, nl_t nl, promise_t update, promise_new_particles_t newParticles) {

        double localAlpha = 0.0;
        std::vector<event_t> localEvents{};
        std::set<event_t> boxoverlappingEvents {};
        gatherEvents<false>(kernel, box.particleIndices, nl, data, localAlpha, approximateRate, localEvents);
        {
            // handle events
            newParticles.set_value(handleEventsGillespie(kernel, false, approximateRate, std::move(localEvents)));
        }
        update.set_value(std::move(boxoverlappingEvents));
    };

    std::vector<std::future<std::set<unsigned long>>> updates;
    std::vector<std::future<std::vector<particle_t>>> newParticles;
    {
        //readdy::util::Timer t ("\t run threads");
        std::vector<util::ScopedThread> threads;
        for (unsigned int i = 0; i < kernel->getNThreads(); ++i) {
            promise_t promise;
            updates.push_back(promise.get_future());
            promise_new_particles_t promiseParticles;
            newParticles.push_back(promiseParticles.get_future());
            threads.push_back(
                    util::ScopedThread(std::thread(
                            worker, std::ref(boxes[i]),
                            std::ref(kernel->getKernelContext()),
                            kernel->getKernelStateModel().getParticleData(),
                            kernel->getKernelStateModel().getNeighborList(),
                            std::move(promise),
                            std::move(promiseParticles)
                    ))
            );
        }
    }
    {
        //readdy::util::Timer t ("\t fix marked");
        std::vector<event_t> evilEvents{};
        double alpha = 0;
        long n_local_problematic = 0;
        for (auto &&update : updates) {
            auto &&local_problematic = update.get();
            n_local_problematic += local_problematic.size();
            gatherEvents<false>(kernel, std::move(local_problematic), kernel->getKernelStateModel().getNeighborList(),
                                kernel->getKernelStateModel().getParticleData(), alpha, approximateRate, evilEvents);
        }
        //BOOST_LOG_TRIVIAL(debug) << "got n_local_problematic="<<n_local_problematic<<", handling events on these!";
        auto newProblemParticles = handleEventsGillespie(kernel, false, approximateRate, std::move(evilEvents));
        // BOOST_LOG_TRIVIAL(trace) << "got problematic particles by conflicts within box: " << n_local_problematic;

        kernel->getKernelStateModel().getParticleData()->deactivateMarked();

        const auto &fixPos = kernel->getKernelContext().getFixPositionFun();
        for (auto &&future : newParticles) {
            auto particles = future.get();
            // reposition particles to respect the periodic b.c.
            std::for_each(particles.begin(), particles.end(),
                          [&fixPos](readdy::model::Particle &p) { fixPos(p.getPos()); });
            kernel->getKernelStateModel().getParticleData()->addParticles(particles);
        }
        std::for_each(newProblemParticles.begin(), newProblemParticles.end(),
                      [&fixPos](readdy::model::Particle &p) { fixPos(p.getPos()); });
        kernel->getKernelStateModel().getParticleData()->addParticles(newProblemParticles);
    }
}
