/**
 * << detailed description >>
 *
 * @file GillespieParallel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#include <future>
#include <queue>

#include <readdy/kernel/cpu/programs/reactions/GillespieParallel.h>


using rdy_particle_t = readdy::model::Particle;

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {


long GillespieParallel::SlicedBox::getShellIndex(const vec_t &pos) const {
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

GillespieParallel::SlicedBox::SlicedBox(unsigned int id, vec_t lowerLeftVertex, vec_t upperRightVertex,
                                        double maxReactionRadius,
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

bool GillespieParallel::SlicedBox::isInBox(const vec_t &particle) const {
    return particle[longestAxis] >= leftBoundary && particle[longestAxis] < rightBoundary;
}

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
    for (auto &e : particleData->entries) {
        if (!e.is_deactivated()) {
            unsigned int boxIndex = static_cast<unsigned int>(
                    floor((e.pos[longestAxis] + .5 * simBoxSize[longestAxis]) / boxWidth)
            );
            if (boxIndex < nBoxes) {
                auto &box = boxes[boxIndex];
                box.particleIndices.push_back(&e);
            }
        }
    }
}

void GillespieParallel::clear() {
    maxReactionRadius = 0;
    boxes.clear();
}

void GillespieParallel::handleBoxReactions() {
    using promise_t = std::promise<std::set<data_t::Entry *>>;
    using promise_new_particles_t = std::promise<data_t::update_t>;

    auto worker = [this](SlicedBox &box, ctx_t ctx, data_t *data, nl_t nl, promise_t update,
                         promise_new_particles_t newParticles) {
        const auto &fixPos = kernel->getKernelContext().getFixPositionFun();
        std::set<data_t::Entry *> problematic{};
        double localAlpha = 0.0;
        std::vector<event_t> localEvents{};
        // step 1: find all problematic particles (ie the ones, that have (also transitively)
        // a reaction with a particle in another box)
        {
            for (const auto pIdx : box.particleIndices) {
                findProblematicParticles(pIdx, box, ctx, *data, nl, problematic);
            }
        }
        // step 2: remove the problematic ones out of the box and gather remaining events
        {
            box.particleIndices.erase(
                    std::remove_if(box.particleIndices.begin(), box.particleIndices.end(),
                                   [&problematic](data_t::Entry *x) {
                                       return problematic.find(x) != problematic.end();
                                   }), box.particleIndices.end()
            );
            gatherEvents(kernel, box.particleIndices, nl, *data, localAlpha, localEvents);
            // handle events
            {
                auto result = handleEventsGillespie(kernel, false, approximateRate, std::move(localEvents));
                std::for_each(std::get<0>(result).begin(), std::get<0>(result).end(),
                              [&fixPos](data_t::Entry &p) { fixPos(p.pos); });
                newParticles.set_value(std::move(result));
            }

        }
        update.set_value(std::move(problematic));
    };

    std::vector<std::future<std::set<data_t::Entry *>>> updates;
    std::vector<std::future<data_t::update_t>> newParticles;
    {
        //readdy::util::Timer t ("\t run threads");
        std::vector<util::scoped_thread> threads;
        for (unsigned int i = 0; i < kernel->getNThreads(); ++i) {
            // nboxes == nthreads
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
            gatherEvents(kernel, std::move(local_problematic), kernel->getKernelStateModel().getNeighborList(),
                         *kernel->getKernelStateModel().getParticleData(), alpha, evilEvents);
        }
        //BOOST_LOG_TRIVIAL(debug) << "got n_local_problematic="<<n_local_problematic<<", handling events on these!";
        auto newProblemParticles = evilEvents.empty() ? data_t::update_t() :
                                   handleEventsGillespie(kernel, false, approximateRate, std::move(evilEvents));
        // BOOST_LOG_TRIVIAL(trace) << "got problematic particles by conflicts within box: " << n_local_problematic;

        auto &data = *kernel->getKernelStateModel().getParticleData();
        auto &neighbor_list = *kernel->getKernelStateModel().getNeighborList();
        const auto &fixPos = kernel->getKernelContext().getFixPositionFun();
        for (auto &&future : newParticles) {
            auto particles = future.get();
            neighbor_list.updateData(data, std::move(particles));
        }
        std::for_each(std::get<0>(newProblemParticles).begin(), std::get<0>(newProblemParticles).end(),
                      [&fixPos](data_t::Entry &p) { fixPos(p.pos); });
        neighbor_list.updateData(data, std::move(newProblemParticles));
    }

}

void GillespieParallel::findProblematicParticles(
        data_t::Entry *entry, const SlicedBox &box, ctx_t ctx,
        const data_t &data, nl_t nl, std::set<data_t::Entry *> &problematic
) const {
    if (problematic.find(entry) != problematic.end()) {
        return;
    }

    // we only want out most particles here, since the transitively dependent particles
    // are resolved within this method and therefore have no significance as input parameter
    if (box.getShellIndex(entry->pos) > 0) {
        //BOOST_LOG_TRIVIAL(debug) << "--> ignoring particle " << (*data)[idx] << " with shell idx " << box.getShellIndex(pPos);
        return;
    }

    std::queue<decltype(entry)> bfs{};

    try {
        for (const auto &neighbor : nl->neighbors(entry)) {
            const auto &reactions = ctx.getOrder2Reactions(entry->type, neighbor.idx->type);
            if (!reactions.empty()) {
                const auto distSquared = neighbor.d2;
                for (const auto &r : reactions) {
                    if (r->getRate() > 0 && distSquared < r->getEductDistanceSquared()) {
                        const bool neighborProblematic = problematic.find(neighbor.idx) != problematic.end();
                        if (neighborProblematic || !box.isInBox(neighbor.idx->pos)) {
                            // we have a problematic particle!
                            problematic.insert(entry);
                            bfs.push(entry);
                            break;
                        }
                    }
                }
            }
        }
    } catch(const std::out_of_range& e) {
        log::console()->trace("The particle {} was not in the neighbor list 2! {}", data.getEntryIndex(entry), e.what());
    }
    //BOOST_LOG_TRIVIAL(debug) << "------------------ BFS -----------------";
    while (!bfs.empty()) {
        const auto x = bfs.front();
        bfs.pop();
        const auto x_shell_idx = box.getShellIndex(x->pos);
        //BOOST_LOG_TRIVIAL(debug) << " ----> looking at neighbors of " << x;
        try {

            for (const auto &x_neighbor : nl->neighbors(x)) {
                //BOOST_LOG_TRIVIAL(debug) << "\t ----> neighbor " << x_neighbor;
                /*if(neighbor_shell_idx == 0) {
                    //BOOST_LOG_TRIVIAL(debug) << "\t ----> neighbor was in outer shell, ignore";
                    continue;
                } else {
                    //BOOST_LOG_TRIVIAL(debug) << "\t ----> got neighbor with shell index " << neighbor_shell_idx;
                }*/
                //BOOST_LOG_TRIVIAL(debug) << "\t\t inBox=" <<box.isInBox(neighborPos) <<", shellabsdiff=" <<std::abs(x_shell_idx - neighbor_shell_idx);
                if (box.isInBox(x_neighbor.idx->pos)
                    && std::abs(x_shell_idx - box.getShellIndex(x_neighbor.idx->pos)) <= 1.0) {
                    //BOOST_LOG_TRIVIAL(debug) << "\t\t neighbor was in box and adjacent shell";
                    const auto &reactions = ctx.getOrder2Reactions(x->type, x_neighbor.idx->type);
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
        } catch (const std::out_of_range& e) {
            log::console()->trace("the particle {} was not in the neighbor list 3! {}", data.getEntryIndex(x), e.what());
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

void GillespieParallel::setApproximateRate(bool approximateRate) {
    GillespieParallel::approximateRate = approximateRate;
}

GillespieParallel::~GillespieParallel() = default;
}

}
}
}
}