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
                    Gillespie::handleEvents(std::vector<Gillespie::_event_t> events, double alpha) {
                        using _rdy_particle_t = readdy::model::Particle;
                        std::vector<_rdy_particle_t> newParticles{};

                        const auto& ctx = kernel->getKernelContext();
                        auto rnd = readdy::model::RandomProvider();
                        auto data = kernel->getKernelStateModel().getParticleData();
                        const auto dt = ctx.getTimeStep();
                        /**
                         * Handle gathered reaction events
                         */
                        {
                            std::size_t nDeactivated = 0;
                            const std::size_t nEvents = events.size();
                            while (nDeactivated < nEvents) {
                                alpha = (*(events.end() - nDeactivated - 1)).cumulativeRate;
                                const auto x = rnd.getUniform(0, alpha);
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
                                if(rnd.getUniform() < event.reactionRate*dt) {
                                    /**
                                     * Perform reaction
                                     */
                                    {
                                        const auto p1 = data->operator[](event.idx1);
                                        _rdy_particle_t pOut1{}, pOut2{};
                                        if (event.nEducts == 1) {
                                            auto reaction = ctx.getOrder1Reactions(event.t1)[event.reactionIdx];
                                            if (reaction->getNProducts() == 1) {
                                                reaction->perform(p1, p1, pOut1, pOut2);
                                                newParticles.push_back(pOut1);
                                            } else if (reaction->getNProducts() == 2) {
                                                reaction->perform(p1, data->operator[](event.idx2), pOut1, pOut2);
                                                newParticles.push_back(pOut1);
                                                newParticles.push_back(pOut2);
                                            }
                                        } else {
                                            auto reaction = ctx.getOrder2Reactions(event.t1,
                                                                                   event.t2)[event.reactionIdx];
                                            const auto p2 = data->operator[](event.idx2);
                                            if (reaction->getNProducts() == 1) {
                                                reaction->perform(p1, p2, pOut1, pOut2);
                                                newParticles.push_back(pOut1);
                                            } else if (reaction->getNProducts() == 2) {
                                                reaction->perform(p1, p2, pOut1, pOut2);
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
                                            data->markForDeactivation((size_t) idx1);
                                            while (_it < events.end() - nDeactivated) {
                                                if ((*_it).idx1 == idx1 ||
                                                    ((*_it).nEducts == 2 && (*_it).idx2 == idx1)) {
                                                    nDeactivated++;
                                                    std::iter_swap(_it, events.end() - nDeactivated);
                                                }
                                                cumsum += (*_it).reactionRate;
                                                (*_it).cumulativeRate = cumsum;
                                                ++_it;
                                            }
                                        } else {
                                            const auto idx2 = event.idx2;
                                            data->markForDeactivation((size_t) event.idx1);
                                            data->markForDeactivation((size_t) event.idx2);
                                            while (_it < events.end() - nDeactivated) {
                                                if ((*_it).idx1 == idx1 || (*_it).idx1 == idx2
                                                    ||
                                                    ((*_it).nEducts == 2 &&
                                                     ((*_it).idx2 == idx1 || (*_it).idx2 == idx2))) {
                                                    nDeactivated++;
                                                    std::iter_swap(_it, events.end() - nDeactivated);
                                                }
                                                cumsum += (*_it).reactionRate;
                                                (*_it).cumulativeRate = cumsum;
                                                ++_it;
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
                        return newParticles;
                    }

                    struct GillespieParallel::HaloBox  {
                        std::vector<HaloBox *> neighboringBoxes{};
                        std::vector<unsigned long> particleIndices{};
                        unsigned int id = 0;
                        vec_t lowerLeftVertex, upperRightVertex;
                        vec_t lowerLeftVertexHalo, upperRightVertexHalo;


                        HaloBox(unsigned int id, vec_t lowerLeftBdry, vec_t upperRightBdry, double haloWidth, unsigned int longestAxis, bool periodicInLongestAx) :
                                id(id), lowerLeftVertexHalo(lowerLeftBdry), upperRightVertexHalo(upperRightBdry) {
                            lowerLeftVertex = lowerLeftVertexHalo;
                            if(periodicInLongestAx && id == 0) lowerLeftVertex[longestAxis] += haloWidth;
                            upperRightVertex = upperRightVertexHalo;
                            if(periodicInLongestAx && id == util::getNThreads()-1) upperRightVertex[longestAxis] -= haloWidth;
                        }

                        void addNeighbor(HaloBox *box) {
                            if (box && box->id != id) neighboringBoxes.push_back(box);
                        }

                        friend bool operator==(const HaloBox &lhs, const HaloBox &rhs) {
                            return lhs.id == rhs.id;
                        }

                        friend bool operator!=(const HaloBox &lhs, const HaloBox &rhs) {
                            return !(lhs == rhs);
                        }

                        bool isInBox(const vec_t& particle) const {
                            return particle >= lowerLeftVertex && particle <= upperRightVertex;
                        }

                        bool isInBox(const unsigned long particleIdx) const {
                            return std::find(particleIndices.begin(), particleIndices.end(), particleIdx) != particleIndices.end();
                        }
                    };

                    long positive_modulo(long i, long n) {
                        return (i % n + n) % n;
                    }

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
                            unsigned int nBoxes = static_cast<unsigned int>(util::getNThreads());
                            auto boxBoundingVertices = kernel->getKernelContext().getBoxBoundingVertices();
                            auto lowerLeft = std::get<0>(boxBoundingVertices);
                            const bool periodicInLongestAx = kernel->getKernelContext().getPeriodicBoundary()[longestAxis];
                            const auto boxWidth = simBoxSize[longestAxis] / util::getNThreads();
                            auto upperRight = lowerLeft;
                            {
                                upperRight[otherAxis1] = std::get<1>(boxBoundingVertices)[otherAxis1];
                                upperRight[otherAxis2] = std::get<1>(boxBoundingVertices)[otherAxis2];
                                upperRight[longestAxis] = lowerLeft[longestAxis] + boxWidth;
                            }
                            boxes.reserve(nBoxes);
                            for(unsigned int i = 0; i < nBoxes; ++i) {
                                HaloBox box {i, lowerLeft, upperRight, .5*maxReactionRadius, longestAxis, periodicInLongestAx};
                                boxes.push_back(std::move(box));
                                lowerLeft[longestAxis] += boxWidth;
                                upperRight[longestAxis] += boxWidth;
                            }

                            for(unsigned int i = 1; i < nBoxes-1; ++i) {
                                boxes[i].addNeighbor(&boxes[positive_modulo(i-1, nBoxes)]);
                                boxes[i].addNeighbor(&boxes[positive_modulo(i+1, nBoxes)]);
                            }
                            if(nBoxes > 1) {
                                boxes[0].addNeighbor(&boxes[1]);
                                boxes[nBoxes-1].addNeighbor(&boxes[nBoxes-2]);
                                if(kernel->getKernelContext().getPeriodicBoundary()[longestAxis]) {
                                    boxes[0].addNeighbor(&boxes[nBoxes-1]);
                                    boxes[nBoxes-1].addNeighbor(&boxes[0]);
                                }
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
                                if(box.isInBox(*it)) {
                                    box.particleIndices.push_back(idx);
                                } else {
                                    problematicParticles.push_back(idx);
                                }
                            }
                            ++idx;
                        }
                    }

                    void GillespieParallel::clear() {
                        maxReactionRadius = 0;
                        problematicParticles.clear();
                        boxes.clear();
                    }

                    void GillespieParallel::handleBoxReactions() {
                        using promise_t = std::promise<std::vector<unsigned long>>;

                        auto worker = [](const HaloBox &box, const ctx_t ctx, data_t data, nl_t nl, promise_t update) {
                            std::vector<unsigned long> problematic;
                            const auto &dist = ctx.getDistSquaredFun();
                            // step 1: find all problematic particles (ie the ones, that have (also transitively)
                            // a reaction with a particle in the halo region
                            {
                                auto it = box.particleIndices.begin();
                                while (it < box.particleIndices.end()) {
                                    handleProblematic(*it, box, ctx, data, nl, problematic);
                                    ++it;
                                }
                            }
                            update.set_value(std::move(problematic));
                        };

                        std::vector<std::future<std::vector<unsigned long>>> updates;
                        {
                            std::vector<util::ScopedThread> threads;
                            for (unsigned int i = 0; i < util::getNThreads(); ++i) {
                                promise_t promise;
                                updates.push_back(promise.get_future());
                                threads.push_back(
                                        util::ScopedThread(std::thread(
                                                worker, std::ref(boxes[i]),
                                                std::ref(kernel->getKernelContext()),
                                                kernel->getKernelStateModel().getParticleData(),
                                                kernel->getKernelStateModel().getNeighborList(),
                                                std::move(promise)
                                        ))
                                );
                            }
                        }
                        for(auto&& update : updates) {
                            BOOST_LOG_TRIVIAL(debug) << "got update: " << update.get()[0];
                        }

                    }

                    void GillespieParallel::handleProblematic(
                            const unsigned long idx, const HaloBox &box, const ctx_t ctx,
                            data_t data, nl_t nl, std::vector<unsigned long> &update
                    ) const {
                        if(std::find(update.begin(), update.end(), idx) != update.end()) return;

                        const auto &dist = ctx.getDistSquaredFun();
                        const auto pType = *(data->begin_types() + idx);
                        const auto pPos = *(data->begin_positions() + idx);
                        auto nlIt = nl->pairs->find(idx);
                        if (nlIt != nl->pairs->end()) {
                            for (const auto neighborIdx : nlIt->second) {
                                const auto neighborType = *(data->begin_types() + neighborIdx);
                                const auto neighborPos = *(data->begin_positions() + neighborIdx);
                                const auto& reactions = ctx.getOrder2Reactions(pType, neighborType);
                                if(!reactions.empty()) {
                                    const auto distSquared = dist(pPos, neighborPos);
                                    for (auto itReactions = reactions.begin(); itReactions != reactions.end(); ++itReactions) {
                                        const auto reaction = *itReactions;
                                        if (distSquared < reaction->getEductDistance() * reaction->getEductDistance()) {
                                            if(reaction->getRate() > 0) {
                                                const bool alreadyProblematic = std::find(update.begin(), update.end(), neighborIdx) != update.end();
                                                if (alreadyProblematic || !box.isInBox(neighborPos)) {
                                                    // we have a problematic particle!
                                                    update.push_back(idx);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    GillespieParallel::~GillespieParallel() = default;
                }

            }
        }
    }
}