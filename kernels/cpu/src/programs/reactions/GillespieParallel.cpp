/**
 * << detailed description >>
 *
 * @file GillespieParallel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#include <readdy/kernel/cpu/programs/reactions/GillespieParallel.h>
#include <future>
#include <readdy/kernel/cpu/programs/reactions/ReactionUtils.h>
#include <queue>


using rdy_particle_t = readdy::model::Particle;

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {


struct GillespieParallel::SlicedBox {
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

    long getShellIndex(const vec_t &pos) const {
        if (shellWidth > 0) {
            const auto mindist = std::min(
                    std::abs(pos[longestAxis] - leftBoundary),
                    std::abs(rightBoundary - pos[longestAxis])
            );
            return static_cast<long>(std::floor(mindist / shellWidth));
        } else {
            return 0;
        }
    }

    SlicedBox(unsigned int id, vec_t lowerLeftVertex, vec_t upperRightVertex, double maxReactionRadius,
              unsigned int longestAxis)
            : id(id), lowerLeftVertex(lowerLeftVertex), upperRightVertex(upperRightVertex), longestAxis(longestAxis) {
        leftBoundary = lowerLeftVertex[longestAxis];
        rightBoundary = upperRightVertex[longestAxis];
        boxWidth = rightBoundary - leftBoundary;
        n_shells = static_cast<particle_indices_t::size_type>(
                std::floor(.5 * boxWidth / maxReactionRadius)
        );
        shellWidth = .5 * boxWidth / static_cast<double>(n_shells);
        particleIndices.resize(n_shells);
    }

    friend bool operator==(const SlicedBox &lhs, const SlicedBox &rhs) {
        return lhs.id == rhs.id;
    }

    friend bool operator!=(const SlicedBox &lhs, const SlicedBox &rhs) {
        return !(lhs == rhs);
    }

    bool isInBox(const vec_t &particle) const {
        return particle[longestAxis] >= leftBoundary && particle[longestAxis] < rightBoundary;
    }
};

void GillespieParallel::execute() {
    {
        //readdy::util::Timer t ("gillespie parallel: setup boxes");
        setupBoxes();
    }
    {
        //readdy::util::Timer t ("gillespie parallel: fill boxes");
        fillBoxes();
    }
    {
        //readdy::util::Timer t ("gillespie parallel: handle reactions");
        handleBoxReactions();
    }
}

GillespieParallel::GillespieParallel(const kernel_t *const kernel) : kernel(kernel), boxes({}) {}

void GillespieParallel::setupBoxes() {
    if (boxes.empty()) {
        double maxReactionRadius = 0.0;
        for (auto &&e : kernel->getKernelContext().getAllOrder2Reactions()) {
            maxReactionRadius = std::max(maxReactionRadius, e->getEductDistance());
        }

        const auto &simBoxSize = kernel->getKernelContext().getBoxSize();
        unsigned int longestAxis{
                static_cast<unsigned int>(
                        std::max_element(simBoxSize.begin(), simBoxSize.end()) - simBoxSize.begin()
                )
        };
        unsigned int otherAxis1 = [longestAxis]() -> unsigned int {
            switch (longestAxis) {
                case 0:
                    return 1;
                case 1:
                    return 2;
                default:
                    return 0;
            }
        }();
        unsigned int otherAxis2 = [longestAxis]() -> unsigned int {
            switch (longestAxis) {
                case 0:
                    return 2;
                case 1:
                    return 0;
                default:
                    return 1;
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
        for (unsigned int i = 0; i < nBoxes; ++i) {
            SlicedBox box{i, lowerLeft, upperRight, maxReactionRadius, longestAxis};
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
    std::for_each(boxes.begin(), boxes.end(), [](SlicedBox &box) { box.particleIndices.clear(); });
    const auto particleData = kernel->getKernelStateModel().getParticleData();
    const auto simBoxSize = kernel->getKernelContext().getBoxSize();
    const auto nBoxes = boxes.size();
    std::size_t idx = 0;
    for (auto it = particleData->begin_positions(); it < particleData->end_positions(); ++it) {
        unsigned int boxIndex = static_cast<unsigned int>(
                floor(((*it)[longestAxis] + .5 * simBoxSize[longestAxis]) / boxWidth)
        );
        if (boxIndex < nBoxes) {
            auto &box = boxes[boxIndex];
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

    auto worker = [this](SlicedBox &box, ctx_t ctx, data_t data, nl_t nl, promise_t update,
                         promise_new_particles_t newParticles) {
        std::set<unsigned long> problematic{};
        double localAlpha = 0.0;
        std::vector<event_t> localEvents{};
        // step 1: find all problematic particles (ie the ones, that have (also transitively)
        // a reaction with a particle in another box)
        {
            for (const auto pIdx : box.particleIndices) {
                findProblematicParticles(pIdx, box, ctx, data, nl, problematic);
            }
        }
        // step 2: remove the problematic ones out of the box and gather remaining events
        {
            box.particleIndices.erase(
                    std::remove_if(box.particleIndices.begin(), box.particleIndices.end(),
                                   [&problematic](const unsigned long x) {
                                       return problematic.find(x) != problematic.end();
                                   }), box.particleIndices.end()
            );
            gatherEvents(box.particleIndices, nl, data, localAlpha, localEvents);
            // handle events
            {
                newParticles.set_value(
                        handleEventsGillespie(kernel, filterEventsInAdvance, approximateRate, std::move(localEvents))
                );
            }

        }
        update.set_value(std::move(problematic));
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
            gatherEvents(std::move(local_problematic), kernel->getKernelStateModel().getNeighborList(),
                         kernel->getKernelStateModel().getParticleData(), alpha, evilEvents);
        }
        //BOOST_LOG_TRIVIAL(debug) << "got n_local_problematic="<<n_local_problematic<<", handling events on these!";
        auto newProblemParticles = handleEventsGillespie(kernel, filterEventsInAdvance, approximateRate,
                                                         std::move(evilEvents));
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

void GillespieParallel::findProblematicParticles(
        const unsigned long idx, const SlicedBox &box, ctx_t ctx,
        data_t data, nl_t nl, std::set<unsigned long> &problematic
) const {
    if (problematic.find(idx) != problematic.end()) {
        return;
    }

    const auto pPos = *(data->begin_positions() + idx);
    // we only want out most particles here, since the transitively dependent particles
    // are resolved within this method and therefore have no significance as input parameter
    if (box.getShellIndex(pPos) > 0) {
        //BOOST_LOG_TRIVIAL(debug) << "--> ignoring particle " << (*data)[idx] << " with shell idx " << box.getShellIndex(pPos);
        return;
    }

    std::queue<unsigned long> bfs{};

    auto nlIt = nl->pairs->find(idx);
    if (nlIt != nl->pairs->end()) {
        const auto pType = *(data->begin_types() + idx);
        for (const auto &neighbor : nlIt->second) {
            const auto neighborType = *(data->begin_types() + neighbor.idx);
            const auto &reactions = ctx.getOrder2Reactions(pType, neighborType);
            if (!reactions.empty()) {
                const auto neighborPos = *(data->begin_positions() + neighbor.idx);
                const auto distSquared = neighbor.d2;
                for (const auto &r : reactions) {
                    if (r->getRate() > 0 && distSquared < r->getEductDistanceSquared()) {
                        const bool neighborProblematic = problematic.find(neighbor.idx) != problematic.end();
                        if (neighborProblematic || !box.isInBox(neighborPos)) {
                            // we have a problematic particle!
                            problematic.insert(idx);
                            bfs.push(idx);
                            break;
                        }
                    }
                }
            }
        }
    }
    //BOOST_LOG_TRIVIAL(debug) << "------------------ BFS -----------------";
    while (!bfs.empty()) {
        const auto x = bfs.front();
        bfs.pop();
        const auto xPos = *(data->begin_positions() + x);
        const auto x_shell_idx = box.getShellIndex(xPos);
        auto neighIt = nl->pairs->find(x);
        if (neighIt != nl->pairs->end()) {
            //BOOST_LOG_TRIVIAL(debug) << " ----> looking at neighbors of " << x;
            const auto x_type = *(data->begin_types() + x);
            for (const auto &x_neighbor : neighIt->second) {
                //BOOST_LOG_TRIVIAL(debug) << "\t ----> neighbor " << x_neighbor;
                const auto neighborType = *(data->begin_types() + x_neighbor.idx);
                const auto neighborPos = *(data->begin_positions() + x_neighbor.idx);
                const auto neighbor_shell_idx = box.getShellIndex(neighborPos);
                /*if(neighbor_shell_idx == 0) {
                    //BOOST_LOG_TRIVIAL(debug) << "\t ----> neighbor was in outer shell, ignore";
                    continue;
                } else {
                    //BOOST_LOG_TRIVIAL(debug) << "\t ----> got neighbor with shell index " << neighbor_shell_idx;
                }*/
                //BOOST_LOG_TRIVIAL(debug) << "\t\t inBox=" <<box.isInBox(neighborPos) <<", shellabsdiff=" <<std::abs(x_shell_idx - neighbor_shell_idx);
                if (box.isInBox(neighborPos) && std::abs(x_shell_idx - neighbor_shell_idx) <= 1.0) {
                    //BOOST_LOG_TRIVIAL(debug) << "\t\t neighbor was in box and adjacent shell";
                    const auto &reactions = ctx.getOrder2Reactions(x_type, neighborType);
                    if (!reactions.empty()) {
                        //BOOST_LOG_TRIVIAL(debug) << "\t\t neighbor had potentially conflicting reactions";
                        const bool alreadyProblematic = problematic.find(x_neighbor.idx) != problematic.end();
                        if (alreadyProblematic) {
                            //BOOST_LOG_TRIVIAL(debug) << "\t\t neighbor was already found, ignore";
                            continue;
                        }
                        const auto distSquared = x_neighbor.d2;
                        for (const auto &reaction : reactions) {
                            if (reaction->getRate() > 0 && distSquared < reaction->getEductDistanceSquared()) {
                                // we have a problematic particle!
                                problematic.insert(x_neighbor.idx);
                                bfs.push(x_neighbor.idx);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

}

template<typename ParticleCollection>
void
GillespieParallel::gatherEvents(const ParticleCollection &particles, const nl_t nl, const data_t data, double &alpha,
                                std::vector<GillespieParallel::event_t> &events) const {
    const auto dt = kernel->getKernelContext().getTimeStep();
    for (const auto idx : particles) {
        // this being false should really not happen, though
        if (!*(data->begin_deactivated() + idx)) {
            const auto particleType = *(data->begin_types() + idx);
            // order 1
            {
                const auto &reactions = kernel->getKernelContext().getOrder1Reactions(particleType);
                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                    const auto rate = (*it)->getRate();
                    if (rate > 0 && (!filterEventsInAdvance || performEvent(rate, dt, approximateRate))) {
                        alpha += rate;
                        events.push_back(
                                {1, (*it)->getNProducts(), idx, 0, rate, alpha,
                                 static_cast<index_t>(it - reactions.begin()),
                                 particleType, 0});
                    }
                }
            }
            // order 2
            {
                auto nl_it = nl->pairs->find(idx);
                if (nl_it != nl->pairs->end()) {
                    for (const auto &idx_neighbor : nl_it->second) {
                        if (idx > idx_neighbor.idx) continue;
                        const auto neighborType = *(data->begin_types() + idx_neighbor.idx);
                        const auto &reactions = kernel->getKernelContext().getOrder2Reactions(
                                particleType, neighborType
                        );
                        if (!reactions.empty()) {
                            const auto distSquared = idx_neighbor.d2;
                            for (auto it = reactions.begin(); it < reactions.end(); ++it) {
                                const auto &react = *it;
                                const auto rate = react->getRate();
                                if (rate > 0
                                    && distSquared < react->getEductDistanceSquared()
                                    && (!filterEventsInAdvance || performEvent(rate, dt, approximateRate))) {
                                    if (*(data->begin_deactivated() + idx_neighbor.idx)) {
                                        auto get_box_idx = [&](
                                                unsigned long _idx) {
                                            return static_cast<unsigned int>( floor(
                                                    ((*(data->begin_positions() + _idx))[longestAxis] +
                                                     .5 * kernel->getKernelContext().getBoxSize()[longestAxis]) /
                                                    boxWidth));
                                        };
                                        const auto nBoxIdx = get_box_idx(idx_neighbor.idx);
                                        const auto nPos = *(data->begin_positions() + idx_neighbor.idx);
                                        const auto nShellIdx = boxes[nBoxIdx].getShellIndex(nPos);
                                        const auto myBoxIdx = get_box_idx(idx);
                                        const auto myPos = *(data->begin_positions() + idx);
                                        log::console()->error("Got neighbor that was deactivated (ie conflicting "
                                                                      "reaction) in shell index {}, was in box[{}]={}",
                                                              nShellIdx, nBoxIdx, boxes[nBoxIdx].isInBox(nPos));
                                        log::console()->error("The neighbor is neighbor of a particle "
                                                                      "in box[{}]={}, shell={}", myBoxIdx,
                                                              boxes[myBoxIdx].isInBox(myPos),
                                                              boxes[myBoxIdx].getShellIndex(myPos));
                                    }
                                    alpha += rate;
                                    events.push_back({2, react->getNProducts(), idx, idx_neighbor.idx,
                                                      rate, alpha,
                                                      static_cast<index_t>(it - reactions.begin()),
                                                      particleType, neighborType});
                                }
                            }
                        }
                    }
                }
            }
        } else {
            log::console()->error("The particles list which was given to gather events contained a particle that "
                                          "was already deactivated. This should not happen!");
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

void GillespieParallel::setFilterEventsInAdvance(bool filterEventsInAdvance) {
    GillespieParallel::filterEventsInAdvance = filterEventsInAdvance;
}

void GillespieParallel::setApproximateRate(bool approximateRate) {
    GillespieParallel::approximateRate = approximateRate;
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