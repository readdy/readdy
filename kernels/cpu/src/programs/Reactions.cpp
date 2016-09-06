/**
 * << detailed description >>
 *
 * @file CPUDefaultReactionProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/programs/Reactions.h>
#include <future>
#include <queue>

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

                    Gillespie::Gillespie(const CPUKernel *const kernel) : kernel(kernel) {}

                    std::vector<singlecpu::programs::reactions::ReactionEvent> Gillespie::gatherEvents(double &alpha) {
                        using index_t = singlecpu::programs::reactions::ReactionEvent::index_type;
                        std::vector<singlecpu::programs::reactions::ReactionEvent> events;
                        const auto& ctx = kernel->getKernelContext();
                        auto data = kernel->getKernelStateModel().getParticleData();
                        const auto &dist = ctx.getDistSquaredFun();
                        /**
                        * Reactions with one educt
                        */
                        {
                            auto it_type = data->begin_types();
                            const auto end = data->end_types();
                            while (it_type != end) {
                                const auto &reactions = ctx.getOrder1Reactions(*it_type);
                                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                                    const auto rate = (*it)->getRate();
                                    if(rate > 0) {
                                        alpha += rate;
                                        events.push_back({1, (index_t) (it_type - data->begin_types()), 0, rate, alpha,
                                                          (index_t) (it - reactions.begin()), *it_type, 0});
                                    }
                                }
                                ++it_type;
                            }
                        }

                        /**
                         * Reactions with two educts
                         */
                        {
                            const auto *neighborList = kernel->getKernelStateModel().getNeighborList();
                            auto typesBegin = data->begin_types();
                            for (auto &&nl_it = neighborList->pairs->begin(); nl_it != neighborList->pairs->end(); ++nl_it) {
                                const index_t idx1 = nl_it->first;
                                auto neighbors = nl_it->second;
                                for (const auto idx2 : neighbors) {
                                    if (idx1 > idx2) continue;
                                    const auto &reactions = ctx.getOrder2Reactions(
                                            *(data->begin_types() + idx1), *(data->begin_types() + idx2)
                                    );

                                    const auto distSquared = dist(
                                            *(data->begin_positions() + idx1), *(data->begin_positions() + idx2)
                                    );

                                    for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                                        // if close enough
                                        const auto reaction = *it;
                                        if (distSquared < reaction->getEductDistance() * reaction->getEductDistance()) {
                                            const auto rate = reaction->getRate();
                                            if (rate > 0) {
                                                alpha += rate;
                                                events.push_back(
                                                        {2, idx1, idx2, rate, alpha,
                                                         (index_t) (it - reactions.begin()),
                                                         *(typesBegin + idx1), *(typesBegin + idx2)});
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        return events;
                    }

                    std::vector<readdy::model::Particle>
                    handleEventsGillespie(CPUKernel const*const kernel, std::vector<readdy::kernel::singlecpu::programs::reactions::ReactionEvent> events) {
                        using _event_t = readdy::kernel::singlecpu::programs::reactions::ReactionEvent;
                        using _rdy_particle_t = readdy::model::Particle;
                        std::vector<_rdy_particle_t> newParticles{};

                        const auto& ctx = kernel->getKernelContext();
                        auto rnd = std::make_unique<readdy::model::RandomProvider>();
                        auto data = kernel->getKernelStateModel().getParticleData();
                        const auto dt = ctx.getTimeStep();
                        /**
                         * Handle gathered reaction events
                         */
                        std::vector<unsigned long> updateDeactivated {};
                        {
                            std::size_t nDeactivated = 0;
                            const std::size_t nEvents = events.size();
                            while (nDeactivated < nEvents) {
                                const auto alpha = (*(events.end() - nDeactivated - 1)).cumulativeRate;
                                const auto x = rnd->getUniform(0, alpha);
                                const auto eventIt = std::lower_bound(
                                        events.begin(), events.end() - nDeactivated, x,
                                        [](const _event_t &elem1, double elem2) {
                                            return elem1.cumulativeRate < elem2;
                                        }
                                );
                                const auto event = *eventIt;
                                if (eventIt == events.end() - nDeactivated) {
                                    throw std::runtime_error("this should not happen (event not found)");
                                }
                                if(rnd->getUniform() < event.reactionRate*dt) {
                                    /**
                                     * Perform reaction
                                     */
                                    {
                                        const auto p1 = data->operator[](event.idx1);
                                        _rdy_particle_t pOut1{}, pOut2{};
                                        if (event.nEducts == 1) {
                                            auto reaction = ctx.getOrder1Reactions(event.t1)[event.reactionIdx];
                                            if (reaction->getNProducts() == 1) {
                                                reaction->perform(p1, p1, pOut1, pOut2, rnd);
                                                newParticles.push_back(pOut1);
                                            } else if (reaction->getNProducts() == 2) {
                                                reaction->perform(p1, data->operator[](event.idx2), pOut1, pOut2, rnd);
                                                newParticles.push_back(pOut1);
                                                newParticles.push_back(pOut2);
                                            }
                                        } else {
                                            auto reaction = ctx.getOrder2Reactions(event.t1,
                                                                                   event.t2)[event.reactionIdx];
                                            const auto p2 = data->operator[](event.idx2);
                                            if (reaction->getNProducts() == 1) {
                                                reaction->perform(p1, p2, pOut1, pOut2, rnd);
                                                newParticles.push_back(pOut1);
                                            } else if (reaction->getNProducts() == 2) {
                                                reaction->perform(p1, p2, pOut1, pOut2, rnd);
                                                newParticles.push_back(pOut1);
                                                newParticles.push_back(pOut2);
                                            }
                                        }
                                    }
                                    /**
                                     * deactivate events whose educts have disappeared (including the just handled one)
                                     */
                                    {
                                        auto _it = events.begin();
                                        double cumsum = 0.0;
                                        const auto idx1 = event.idx1;
                                        if (event.nEducts == 1) {
                                            *(data->begin_deactivated() + idx1) = true;
                                            updateDeactivated.push_back(idx1);
                                            while (_it < events.end() - nDeactivated) {
                                                if ((*_it).idx1 == idx1 ||
                                                    ((*_it).nEducts == 2 && (*_it).idx2 == idx1)) {
                                                    nDeactivated++;
                                                    std::iter_swap(_it, events.end() - nDeactivated);
                                                    cumsum += (*_it).reactionRate;
                                                    (*_it).cumulativeRate = cumsum;
                                                } else {
                                                    ++_it;
                                                }
                                            }
                                        } else {
                                            const auto idx2 = event.idx2;
                                            *(data->begin_deactivated() + event.idx1) = true;
                                            *(data->begin_deactivated() + event.idx2) = true;
                                            updateDeactivated.push_back(event.idx1);
                                            updateDeactivated.push_back(event.idx2);
                                            while (_it < events.end() - nDeactivated) {
                                                if ((*_it).idx1 == idx1 || (*_it).idx1 == idx2 ||
                                                        ((*_it).nEducts == 2 &&
                                                                ((*_it).idx2 == idx1 || (*_it).idx2 == idx2))) {
                                                    nDeactivated++;
                                                    std::iter_swap(_it, events.end() - nDeactivated);
                                                    cumsum += (*_it).reactionRate;
                                                    (*_it).cumulativeRate = cumsum;
                                                } else {
                                                    ++_it;
                                                }
                                            }
                                        }

                                    }
                                } else {
                                    nDeactivated++;
                                    std::iter_swap(eventIt, events.end() - nDeactivated);
                                    (*eventIt).cumulativeRate = (*eventIt).reactionRate;
                                    if(eventIt > events.begin()) {
                                        (*eventIt).cumulativeRate += (*(eventIt-1)).cumulativeRate;
                                    }
                                    auto cumsum = (*eventIt).cumulativeRate;
                                    for(auto _it = eventIt+1; _it < events.end() - nDeactivated; ++_it) {
                                        cumsum += (*_it).reactionRate;
                                        (*_it).cumulativeRate = cumsum;
                                    }
                                }
                            }
                        }
                        data->updateDeactivated(updateDeactivated);
                        return newParticles;
                    }

                    struct GillespieParallel::HaloBox  {
                        // shellIdx -> particles, outmost shell has idx 0
                        using particle_indices_t = std::vector<unsigned long>;
                        particle_indices_t particleIndices{};
                        unsigned int id = 0;
                        vec_t lowerLeftVertex, upperRightVertex;
                        double leftBoundary = 0;
                        double rightBoundary = 0;
                        particle_indices_t::size_type n_shells;
                        unsigned int longestAxis;
                        double boxWidth;
                        double shellWidth = 0.0;

                        long getShellIndex(const vec_t& pos) const {
                            const auto mindist = std::min(
                                    std::abs(pos[longestAxis] - leftBoundary),
                                    std::abs(pos[longestAxis] - rightBoundary)
                            );
                            return static_cast<long>(std::floor(mindist/shellWidth));
                        }

                        HaloBox(unsigned int id, vec_t lowerLeftVertex, vec_t upperRightVertex, double maxReactionRadius,
                                unsigned int longestAxis)
                                : id(id), lowerLeftVertex(lowerLeftVertex), upperRightVertex(upperRightVertex), longestAxis(longestAxis) {
                            leftBoundary = lowerLeftVertex[longestAxis];
                            rightBoundary = upperRightVertex[longestAxis];
                            boxWidth = rightBoundary - leftBoundary;
                            n_shells = static_cast<particle_indices_t::size_type>(
                                    std::floor(.5*boxWidth/maxReactionRadius)
                            );
                            shellWidth = .5*boxWidth / static_cast<double>(n_shells);
                            particleIndices.resize(n_shells);
                        }

                        friend bool operator==(const HaloBox &lhs, const HaloBox &rhs) {
                            return lhs.id == rhs.id;
                        }

                        friend bool operator!=(const HaloBox &lhs, const HaloBox &rhs) {
                            return !(lhs == rhs);
                        }

                        bool isInBox(const vec_t& particle) const {
                            return particle[longestAxis] >= leftBoundary && particle[longestAxis] < rightBoundary;
                        }
                    };

                    void GillespieParallel::execute() {
                        setupBoxes();
                        fillBoxes();
                        handleBoxReactions();
                    }

                    GillespieParallel::GillespieParallel(const kernel_t *const kernel) : kernel(kernel), boxes({}) { }

                    void GillespieParallel::setupBoxes() {
                        if(boxes.empty()) {
                            double maxReactionRadius = 0.0;
                            for(auto&& e : kernel->getKernelContext().getAllOrder2Reactions()) {
                                maxReactionRadius = std::max(maxReactionRadius, e->getEductDistance());
                            }

                            const auto& simBoxSize = kernel->getKernelContext().getBoxSize();
                            unsigned int longestAxis {
                                    static_cast<unsigned int>(
                                            std::max_element(simBoxSize.begin(), simBoxSize.end()) - simBoxSize.begin()
                                    )
                            };
                            unsigned int otherAxis1 = [longestAxis] () -> unsigned int {
                               switch(longestAxis) {
                                   case 0: return 1;
                                   case 1: return 2;
                                   default: return 0;
                               }
                            }();
                            unsigned int otherAxis2 = [longestAxis] () ->unsigned int {
                                switch(longestAxis) {
                                    case 0: return 2;
                                    case 1: return 0;
                                    default: return 1;
                                }
                            }();
                            unsigned int nBoxes = static_cast<unsigned int>(kernel->getNThreads());
                            auto boxBoundingVertices = kernel->getKernelContext().getBoxBoundingVertices();
                            auto lowerLeft = std::get<0>(boxBoundingVertices);
                            const auto boxWidth = simBoxSize[longestAxis] / nBoxes;
                            auto upperRight = lowerLeft;
                            {
                                upperRight[otherAxis1] = std::get<1>(boxBoundingVertices)[otherAxis1];
                                upperRight[otherAxis2] = std::get<1>(boxBoundingVertices)[otherAxis2];
                                upperRight[longestAxis] = lowerLeft[longestAxis] + boxWidth;
                            }
                            boxes.reserve(nBoxes);
                            for(unsigned int i = 0; i < nBoxes; ++i) {
                                HaloBox box {i, lowerLeft, upperRight, maxReactionRadius, longestAxis};
                                boxes.push_back(std::move(box));
                                lowerLeft[longestAxis] += boxWidth;
                                upperRight[longestAxis] += boxWidth;
                            }

                            GillespieParallel::maxReactionRadius = maxReactionRadius;
                            GillespieParallel::longestAxis = longestAxis;
                            GillespieParallel::otherAxis1 = otherAxis1;
                            GillespieParallel::otherAxis2 = otherAxis2;
                            GillespieParallel::boxWidth = boxWidth;
                        }
                    }

                    void GillespieParallel::fillBoxes() {
                        std::for_each(boxes.begin(), boxes.end(), [](HaloBox& box) {box.particleIndices.clear();});
                        const auto particleData = kernel->getKernelStateModel().getParticleData();
                        const auto simBoxSize = kernel->getKernelContext().getBoxSize();
                        const auto nBoxes = boxes.size();
                        std::size_t idx = 0;
                        for(auto it = particleData->begin_positions(); it < particleData->end_positions(); ++it) {
                            unsigned int boxIndex = static_cast<unsigned int>(
                                    floor(((*it)[longestAxis] + .5 * simBoxSize[longestAxis])/boxWidth)
                            );
                            if(boxIndex < nBoxes) {
                                auto& box = boxes[boxIndex];
                                box.particleIndices.push_back(idx);
                            }
                            ++idx;
                        }
                    }

                    void GillespieParallel::clear() {
                        maxReactionRadius = 0;
                        boxes.clear();
                    }

                    void GillespieParallel::handleBoxReactions() {
                        using promise_t = std::promise<std::set<unsigned long>>;
                        using promise_new_particles_t = std::promise<std::vector<particle_t>>;

                        auto worker = [this](HaloBox &box, ctx_t ctx, data_t data, nl_t nl, promise_t update, promise_new_particles_t newParticles) {
                            std::set<unsigned long> problematic {};
                            double localAlpha = 0.0;
                            std::vector<event_t> localEvents {};
                            // step 1: find all problematic particles (ie the ones, that have (also transitively)
                            // a reaction with a particle in another box)
                            {
                                auto it = box.particleIndices.begin();
                                while (it < box.particleIndices.end()) {
                                    handleProblematic(*it, box, ctx, data, nl, problematic);
                                    ++it;
                                }
                            }
                            // step 2: remove the problematic ones out of the box and gather remaining events
                            {
                                box.particleIndices.erase(
                                        std::remove_if(box.particleIndices.begin(), box.particleIndices.end(), [&problematic](const unsigned long x) {
                                            return problematic.find(x) != problematic.end();
                                        }), box.particleIndices.end()
                                );
                                gatherEvents(box.particleIndices, nl, data, localAlpha, localEvents);
                                // handle events
                                {
                                    newParticles.set_value(
                                            handleEventsGillespie(kernel, std::move(localEvents))
                                    );
                                }

                            }
                            update.set_value(std::move(problematic));
                        };

                        std::vector<std::future<std::set<unsigned long>>> updates;
                        std::vector<std::future<std::vector<particle_t>>> newParticles;
                        {
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
                        std::vector<event_t> evilEvents {};
                        double alpha = 0;
                        long n_local_problematic = 0;
                        for(auto&& update : updates) {
                            auto&& local_problematic = update.get();
                            n_local_problematic += local_problematic.size();
                            gatherEvents(std::move(local_problematic), kernel->getKernelStateModel().getNeighborList(), kernel->getKernelStateModel().getParticleData(), alpha, evilEvents);
                        }
                        // BOOST_LOG_TRIVIAL(trace) << "got problematic particles by conflicts within box: " << n_local_problematic;

                        kernel->getKernelStateModel().getParticleData()->deactivateMarked();

                        const auto &fixPos = kernel->getKernelContext().getFixPositionFun();
                        for(auto&& future : newParticles) {
                            auto particles = future.get();
                            // reposition particles to respect the periodic b.c.
                            std::for_each(particles.begin(), particles.end(),
                                          [&fixPos](readdy::model::Particle &p) { fixPos(p.getPos()); });
                            kernel->getKernelStateModel().getParticleData()->addParticles(particles);
                        }

                    }

                    void GillespieParallel::handleProblematic(
                            const unsigned long idx, const HaloBox &box, ctx_t ctx,
                            data_t data, nl_t nl, std::set<unsigned long> &update
                    ) const {
                        if (update.find(idx) != update.end()) {
                            return;
                        }
                        const auto pPos = *(data->begin_positions() + idx);
                        // we only want out most particles here, since the transitively dependent particles
                        // are resolved within this method and therefore have no significance as input parameter
                        if(box.getShellIndex(pPos) > 0) {
                            //BOOST_LOG_TRIVIAL(debug) << "--> ignoring particle " << (*data)[idx] << " with shell idx " << box.getShellIndex(pPos);
                            return;
                        }

                        std::queue<unsigned long> bfs{};

                        const auto &dist = ctx.getDistSquaredFun();

                        auto nlIt = nl->pairs->find(idx);
                        if (nlIt != nl->pairs->end()) {
                            const auto pType = *(data->begin_types() + idx);
                            for (const auto neighborIdx : nlIt->second) {
                                const auto neighborType = *(data->begin_types() + neighborIdx);
                                const auto neighborPos = *(data->begin_positions() + neighborIdx);
                                const auto &reactions = ctx.getOrder2Reactions(pType, neighborType);
                                if (!reactions.empty()) {
                                    const auto distSquared = dist(pPos, neighborPos);
                                    for (const auto& r : reactions) {
                                        if (r->getRate() > 0) {
                                            if (distSquared < r->getEductDistance() * r->getEductDistance()) {
                                                const bool neighborProblematic = update.find(neighborIdx) != update.end();
                                                if (neighborProblematic || !box.isInBox(neighborPos)) {
                                                    //BOOST_LOG_TRIVIAL(debug) << " ----> found (outer) problematic particle " << (*data)[idx];
                                                    // we have a problematic particle!
                                                    update.insert(idx);
                                                    bfs.push(idx);
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        //BOOST_LOG_TRIVIAL(debug) << "------------------ BFS -----------------";
                        while (!bfs.empty()) {
                            const auto x = bfs.front();
                            const auto xPos = *(data->begin_positions() + x);
                            const auto x_shell_idx = box.getShellIndex(xPos);
                            bfs.pop();
                            auto neighIt = nl->pairs->find(x);
                            if (neighIt != nl->pairs->end()) {
                                //BOOST_LOG_TRIVIAL(debug) << " ----> looking at neighbors of " << x;
                                const auto x_type = *(data->begin_types() + idx);
                                for (const auto x_neighbor : neighIt->second) {
                                    //BOOST_LOG_TRIVIAL(debug) << "\t ----> neighbor " << x_neighbor;
                                    const auto neighborType = *(data->begin_types() + x_neighbor);
                                    const auto neighborPos = *(data->begin_positions() + x_neighbor);
                                    const auto neighbor_shell_idx = box.getShellIndex(neighborPos);
                                    if(neighbor_shell_idx == 0) {
                                        //BOOST_LOG_TRIVIAL(debug) << "\t ----> neighbor was in outer shell, ignore";
                                        continue;
                                    } else {
                                        //BOOST_LOG_TRIVIAL(debug) << "\t ----> got neighbor with shell index " << neighbor_shell_idx;
                                    }
                                    //BOOST_LOG_TRIVIAL(debug) << "\t\t inBox=" <<box.isInBox(neighborPos) <<", shellabsdiff=" <<std::abs(x_shell_idx - neighbor_shell_idx);
                                    if(box.isInBox(neighborPos) && std::abs(x_shell_idx - neighbor_shell_idx) <= 1) {
                                        //BOOST_LOG_TRIVIAL(debug) << "\t\t neighbor was in box and adjacent shell";
                                        const auto &reactions = ctx.getOrder2Reactions(x_type, neighborType);
                                        if (!reactions.empty()) {
                                            //BOOST_LOG_TRIVIAL(debug) << "\t\t neighbor had potentially conflicting reactions";
                                            const bool alreadyProblematic = update.find(x_neighbor) != update.end();
                                            if (alreadyProblematic){
                                                //BOOST_LOG_TRIVIAL(debug) << "\t\t neighbor was already found, ignore";
                                                continue;
                                            }
                                            const auto distSquared = dist(xPos, neighborPos);
                                            for (const auto& reaction : reactions) {
                                                if (distSquared < reaction->getEductDistance() * reaction->getEductDistance()) {
                                                    if (reaction->getRate() > 0) {
                                                        // we have a problematic particle!
                                                        //BOOST_LOG_TRIVIAL(debug) << "\t ----> added neighbor " << x_neighbor<< " to problematic particles";
                                                        update.insert(x_neighbor);
                                                        bfs.push(x_neighbor);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }

                    template<typename ParticleCollection>
                    void GillespieParallel::gatherEvents(const ParticleCollection &particles, const nl_t nl, const data_t data, double &alpha,
                                                         std::vector<GillespieParallel::event_t> &events) const {
                        const auto &dist = kernel->getKernelContext().getDistSquaredFun();
                        for(const auto idx : particles) {
                            const auto particleType = *(data->begin_types() + idx);
                            // order 1
                            {
                                const auto &reactions = kernel->getKernelContext().getOrder1Reactions(particleType);
                                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                                    const auto rate = (*it)->getRate();
                                    if (rate > 0) {
                                        alpha += rate;
                                        events.push_back(
                                                {1, idx, 0, rate, alpha, static_cast<index_t>(it - reactions.begin()),
                                                 particleType, 0});
                                    }
                                }
                            }
                            // order 2
                            {
                                auto nl_it = nl->pairs->find(idx);
                                if(nl_it != nl->pairs->end()) {
                                    for(const auto idx_neighbor : nl_it->second) {
                                        if(idx > idx_neighbor) continue;
                                        const auto neighborType = *(data->begin_types() + idx_neighbor);
                                        const auto& reactions = kernel->getKernelContext().getOrder2Reactions(particleType, neighborType);
                                        if(!reactions.empty()) {
                                            const auto distSquared = dist(
                                                    *(data->begin_positions() + idx), *(data->begin_positions() + idx_neighbor)
                                            );
                                            for(auto it = reactions.begin(); it < reactions.end(); ++it) {
                                                const auto& react = *it;
                                                if(distSquared < react->getEductDistance() * react->getEductDistance()) {
                                                    const auto rate = react->getRate();
                                                    if(rate > 0) {
                                                        alpha += rate;
                                                        events.push_back({2, idx, idx_neighbor, rate, alpha,
                                                                               static_cast<index_t>(it - reactions.begin()),
                                                                               particleType, neighborType});
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }

                    double GillespieParallel::getMaxReactionRadius() const {
                        return maxReactionRadius;
                    }

                    double GillespieParallel::getBoxWidth() const {
                        return boxWidth;
                    }

                    unsigned int GillespieParallel::getLongestAxis() const {
                        return longestAxis;
                    }

                    unsigned int GillespieParallel::getOtherAxis1() const {
                        return otherAxis1;
                    }

                    unsigned int GillespieParallel::getOtherAxis2() const {
                        return otherAxis2;
                    }

                    template void GillespieParallel::gatherEvents<std::vector<unsigned long>>(
                            const std::vector<unsigned long> &particles, const nl_t nl, const data_t data, double &alpha,
                            std::vector<GillespieParallel::event_t> &events) const;

                    template void GillespieParallel::gatherEvents<std::set<unsigned long>>(
                            const std::set<unsigned long> &particles, const nl_t nl, const data_t data, double &alpha,
                            std::vector<GillespieParallel::event_t> &events) const;

                    GillespieParallel::~GillespieParallel() = default;
                }

            }
        }
    }
}