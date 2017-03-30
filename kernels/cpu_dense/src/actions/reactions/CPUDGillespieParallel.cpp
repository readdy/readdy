/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file GillespieParallel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <future>
#include <queue>

#include <readdy/common/thread/scoped_thread.h>

#include <readdy/kernel/cpu_dense/actions/reactions/CPUDGillespieParallel.h>

using particle_t = readdy::model::Particle;

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace actions {
namespace reactions {

namespace thd = readdy::util::thread;


long CPUDGillespieParallel::SlicedBox::getShellIndex(const vec_t &pos) const {
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

CPUDGillespieParallel::SlicedBox::SlicedBox(unsigned int id, vec_t lowerLeftVertex, vec_t upperRightVertex,
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

bool CPUDGillespieParallel::SlicedBox::isInBox(const vec_t &particle) const {
    return particle[longestAxis] >= leftBoundary && particle[longestAxis] < rightBoundary;
}

void CPUDGillespieParallel::perform() {
    {
        setupBoxes();
    }
    {
        fillBoxes();
    }
    {
        handleBoxReactions();
    }
}

CPUDGillespieParallel::CPUDGillespieParallel(kernel_t *const kernel, double timeStep)
        : super(timeStep), kernel(kernel), boxes({}) {}

void CPUDGillespieParallel::setupBoxes() {
    if (boxes.empty()) {
        double maxReactionRadius = 0.0;
        for (auto &&e : kernel->getKernelContext().reactions().order2_flat()) {
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

        CPUDGillespieParallel::maxReactionRadius = maxReactionRadius;
        CPUDGillespieParallel::longestAxis = longestAxis;
        CPUDGillespieParallel::otherAxis1 = otherAxis1;
        CPUDGillespieParallel::otherAxis2 = otherAxis2;
        CPUDGillespieParallel::boxWidth = boxWidth;
    }
}

void CPUDGillespieParallel::fillBoxes() {
    std::for_each(boxes.begin(), boxes.end(), [](SlicedBox &box) { box.particleIndices.clear(); });
    const auto particleData = kernel->getCPUDKernelStateModel().getParticleData();
    const auto simBoxSize = kernel->getKernelContext().getBoxSize();
    const auto nBoxes = boxes.size();
    std::size_t idx = 0;
    for (const auto &e : *particleData) {
        unsigned int boxIndex = static_cast<unsigned int>(
                floor((e.position()[longestAxis] + .5 * simBoxSize[longestAxis]) / boxWidth)
        );
        if (boxIndex < nBoxes) {
            auto &box = boxes[boxIndex];
            box.particleIndices.push_back(idx);
        }
        ++idx;
    }
}

void CPUDGillespieParallel::clear() {
    maxReactionRadius = 0;
    boxes.clear();
}

void CPUDGillespieParallel::handleBoxReactions() {
    using promise_t = std::promise<std::set<data_t::index_t>>;
    using promise_new_particles_t = std::promise<data_t::update_t>;

    auto worker = [this](SlicedBox &box, ctx_t ctx, data_t *data, nl_t nl, promise_t update,
                         promise_new_particles_t newParticles) {
        const auto &fixPos = kernel->getKernelContext().getFixPositionFun();
        const auto &d2 = kernel->getKernelContext().getDistSquaredFun();
        std::set<data_t::index_t> problematic{};
        double localAlpha = 0.0;
        std::vector<event_t> localEvents{};
        // step 1: find all problematic particles (ie the ones, that have (also transitively)
        // a reaction with a particle in another box)
        {
            for (const auto pIdx : box.particleIndices) {
                findProblematicParticles(pIdx, box, ctx, *data, nl, problematic, d2);
            }
        }
        // step 2: remove the problematic ones out of the box and gather remaining events
        {
            box.particleIndices.erase(
                    std::remove_if(box.particleIndices.begin(), box.particleIndices.end(),
                                   [&problematic](data_t::index_t x) {
                                       return problematic.find(x) != problematic.end();
                                   }), box.particleIndices.end()
            );
            gatherEvents(kernel, box.particleIndices, nl, *data, localAlpha, localEvents, d2);
            // handle events
            {
                auto result = handleEventsGillespie(kernel, timeStep, false, approximateRate, std::move(localEvents));
                newParticles.set_value(std::move(result));
            }

        }
        update.set_value(std::move(problematic));
    };

    std::vector<std::future<std::set<data_t::index_t>>> updates;
    std::vector<std::future<data_t::update_t>> newParticles;
    {
        //readdy::util::Timer t ("\t run threads");
        std::vector<thd::scoped_thread> threads;
        for (unsigned int i = 0; i < kernel->getNThreads(); ++i) {
            // nboxes == nthreads
            promise_t promise;
            updates.push_back(promise.get_future());
            promise_new_particles_t promiseParticles;
            newParticles.push_back(promiseParticles.get_future());
            threads.push_back(
                    thd::scoped_thread(std::thread(
                            worker, std::ref(boxes[i]),
                            std::ref(kernel->getKernelContext()),
                            kernel->getCPUDKernelStateModel().getParticleData(),
                            kernel->getCPUDKernelStateModel().getNeighborList(),
                            std::move(promise),
                            std::move(promiseParticles)
                    ))
            );
        }
    }
    {
        //readdy::util::Timer t ("\t fix marked");
        auto &data = *kernel->getCPUDKernelStateModel().getParticleData();
        const auto& d2 = kernel->getKernelContext().getDistSquaredFun();
        auto &neighbor_list = *kernel->getCPUDKernelStateModel().getNeighborList();
        std::vector<event_t> evilEvents{};
        double alpha = 0;
        long n_local_problematic = 0;
        for (auto&& update : updates) {
            auto local_problematic = std::move(update.get());
            n_local_problematic += local_problematic.size();
            gatherEvents(kernel, std::move(local_problematic), &neighbor_list, data, alpha, evilEvents, d2);
        }
        //BOOST_LOG_TRIVIAL(debug) << "got n_local_problematic="<<n_local_problematic<<", handling events on these!";
        auto newProblemParticles = handleEventsGillespie(kernel, timeStep, false, approximateRate, std::move(evilEvents));
        // BOOST_LOG_TRIVIAL(trace) << "got problematic particles by conflicts within box: " << n_local_problematic;

        data.deactivateMarked();
        for (auto &&future : newParticles) {
            data.update(std::move(future.get()));
        }
        data.update(std::move(newProblemParticles));
    }

}

void CPUDGillespieParallel::findProblematicParticles(
        data_t::index_t index, const SlicedBox &box, ctx_t ctx,
        const data_t &data, nl_t nl, std::set<data_t::index_t> &problematic,
        const readdy::model::KernelContext::dist_squared_fun& d2
) const {
    if (problematic.find(index) != problematic.end()) {
        return;
    }

    const auto& me = data.entry_at(index);

    // we only want out most particles here, since the transitively dependent particles
    // are resolved within this method and therefore have no significance as input parameter
    if (box.getShellIndex(me.position()) > 0) {
        //BOOST_LOG_TRIVIAL(debug) << "--> ignoring particle " << (*data)[idx] << " with shell idx " << box.getShellIndex(pPos);
        return;
    }

    std::queue<decltype(index)> bfs{};

    for (const auto& neighbor : nl->find_neighbors(index)) {
        const auto& neighborEntry = data.entry_at(neighbor.idx);
        const auto &reactions = ctx.reactions().order2_by_type(me.type, neighborEntry.type);
        if (!reactions.empty()) {
            const auto distSquared = neighbor.d2;
            for (const auto &r : reactions) {
                if (r->getRate() > 0 && distSquared < r->getEductDistanceSquared()) {
                    const bool neighborProblematic = problematic.find(neighbor.idx) != problematic.end();
                    if (neighborProblematic || !box.isInBox(neighborEntry.position())) {
                        // we have a problematic particle!
                        problematic.insert(index);
                        bfs.push(index);
                        break;
                    }
                }
            }
        }
    }
    //BOOST_LOG_TRIVIAL(debug) << "------------------ BFS -----------------";
    while (!bfs.empty()) {
        const auto x = bfs.front();
        const auto& x_entry = data.entry_at(x);
        bfs.pop();
        const auto x_shell_idx = box.getShellIndex(x_entry.position());
        //BOOST_LOG_TRIVIAL(debug) << " ----> looking at neighbors of " << x;
        for (const auto& x_neighbor : nl->find_neighbors(x)) {
            const auto& x_neighbor_entry = data.entry_at(x_neighbor.idx);
            //BOOST_LOG_TRIVIAL(debug) << "\t ----> neighbor " << x_neighbor;
            /*if(neighbor_shell_idx == 0) {
                //BOOST_LOG_TRIVIAL(debug) << "\t ----> neighbor was in outer shell, ignore";
                continue;
            } else {
                //BOOST_LOG_TRIVIAL(debug) << "\t ----> got neighbor with shell index " << neighbor_shell_idx;
            }*/
            //BOOST_LOG_TRIVIAL(debug) << "\t\t inBox=" <<box.isInBox(neighborPos) <<", shellabsdiff=" <<std::abs(x_shell_idx - neighbor_shell_idx);
            if (box.isInBox(x_neighbor_entry.position())
                && std::abs(x_shell_idx - box.getShellIndex(x_neighbor_entry.position())) <= 1.0) {
                //BOOST_LOG_TRIVIAL(debug) << "\t\t neighbor was in box and adjacent shell";
                const auto &reactions = ctx.reactions().order2_by_type(x_entry.type, x_neighbor_entry.type);
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

double CPUDGillespieParallel::getMaxReactionRadius() const {
    return maxReactionRadius;
}

double CPUDGillespieParallel::getBoxWidth() const {
    return boxWidth;
}

unsigned int CPUDGillespieParallel::getLongestAxis() const {
    return longestAxis;
}

unsigned int CPUDGillespieParallel::getOtherAxis1() const {
    return otherAxis1;
}

unsigned int CPUDGillespieParallel::getOtherAxis2() const {
    return otherAxis2;
}

void CPUDGillespieParallel::setApproximateRate(bool approximateRate) {
    CPUDGillespieParallel::approximateRate = approximateRate;
}

CPUDGillespieParallel::~CPUDGillespieParallel() = default;
}

}
}
}
}