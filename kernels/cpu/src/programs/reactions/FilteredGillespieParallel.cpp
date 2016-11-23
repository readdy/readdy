#include <future>
#include <readdy/kernel/cpu/programs/reactions/FilteredGillespieParallel.h>

/**
 * << detailed description >>
 *
 * @file FilteredGillespieParallel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

using data_t = readdy::kernel::cpu::model::ParticleData;

readdy::kernel::cpu::programs::reactions::FilteredGillespieParallel::FilteredGillespieParallel(
        const readdy::kernel::cpu::Kernel *const kernel) : GillespieParallel(kernel) {}

void readdy::kernel::cpu::programs::reactions::FilteredGillespieParallel::handleBoxReactions() {
    using promise_t = std::promise<std::set<event_t>>;
    using promise_new_particles_t = std::promise<data_t::entries_t>;
    /*

    auto worker = [this](SlicedBox &box, ctx_t ctx, data_t* data, nl_t nl, promise_t update, promise_new_particles_t newParticles) {

        double localAlpha = 0.0;
        std::vector<event_t> localEvents{};
        std::set<event_t> boxoverlappingEvents {};
        gatherEvents(kernel, box.particleIndices, nl, *data, localAlpha, localEvents);
        {
            // handle events
            newParticles.set_value(handleEventsGillespie(kernel, false, approximateRate, std::move(localEvents)));
        }
        update.set_value(std::move(boxoverlappingEvents));
    };

    std::vector<std::future<std::set<event_t>>> updates;
    std::vector<std::future<data_t::entries_t>> newParticles;
    {
        //readdy::util::Timer t ("\t run threads");
        std::vector<util::scoped_thread> threads;
        for (unsigned int i = 0; i < kernel->getNThreads(); ++i) {
            promise_t promise;
            updates.push_back(promise.get_future());
            promise_new_particles_t promiseParticles;
            newParticles.push_back(promiseParticles.get_future());
            threads.push_back(
                    util::scoped_thread(std::thread(
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
        }
        //BOOST_LOG_TRIVIAL(debug) << "got n_local_problematic="<<n_local_problematic<<", handling events on these!";
        auto newProblemParticles = handleEventsGillespie(kernel, false, approximateRate, std::move(evilEvents));
        // BOOST_LOG_TRIVIAL(trace) << "got problematic particles by conflicts within box: " << n_local_problematic;

        const auto &fixPos = kernel->getKernelContext().getFixPositionFun();
        for (auto &&future : newParticles) {
            auto particles = future.get();
            // reposition particles to respect the periodic b.c.
            std::for_each(particles.begin(), particles.end(),
                          [&fixPos](data_t::Entry &p) { fixPos(p.pos); });
            kernel->getKernelStateModel().getParticleData()->addEntries(particles);
        }
        std::for_each(newProblemParticles.begin(), newProblemParticles.end(),
                      [&fixPos](data_t::Entry &p) { fixPos(p.pos); });
        kernel->getKernelStateModel().getParticleData()->addEntries(newProblemParticles);
    }*/
}
