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
 * @date 20.10.16
 */

#include <future>
#include <queue>

#include <readdy/kernel/cpu/actions/reactions/CPUGillespieParallel.h>
#include <readdy/common/thread/scoped_thread.h>


using rdy_particle_t = readdy::model::Particle;

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

namespace thd = readdy::util::thread;


long CPUGillespieParallel::SlicedBox::getShellIndex(const vec_t &pos) const {
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

CPUGillespieParallel::SlicedBox::SlicedBox(unsigned int id, vec_t lowerLeftVertex, vec_t upperRightVertex,
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

bool CPUGillespieParallel::SlicedBox::isInBox(const vec_t &particle) const {
    return particle[longestAxis] >= leftBoundary && particle[longestAxis] < rightBoundary;
}

void CPUGillespieParallel::perform() {
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

CPUGillespieParallel::CPUGillespieParallel(kernel_t *const kernel, double timeStep)
        : super(timeStep), kernel(kernel), boxes({}) {}

void CPUGillespieParallel::setupBoxes() {
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

        CPUGillespieParallel::maxReactionRadius = maxReactionRadius;
        CPUGillespieParallel::longestAxis = longestAxis;
        CPUGillespieParallel::otherAxis1 = otherAxis1;
        CPUGillespieParallel::otherAxis2 = otherAxis2;
        CPUGillespieParallel::boxWidth = boxWidth;
    }
}

void CPUGillespieParallel::fillBoxes() {
    std::for_each(boxes.begin(), boxes.end(), [](SlicedBox &box) { box.particleIndices.clear(); });
    const auto particleData = kernel->getCPUKernelStateModel().getParticleData();
    const auto simBoxSize = kernel->getKernelContext().getBoxSize();
    const auto nBoxes = boxes.size();
    std::size_t idx = 0;
    for (const auto &e : *particleData) {
        if (!e.is_deactivated()) {
            unsigned int boxIndex = static_cast<unsigned int>(
                    floor((e.position()[longestAxis] + .5 * simBoxSize[longestAxis]) / boxWidth)
            );
            if (boxIndex < nBoxes) {
                auto &box = boxes[boxIndex];
                box.particleIndices.push_back(idx);
            }
        }
        ++idx;
    }
}

void CPUGillespieParallel::clear() {
    maxReactionRadius = 0;
    boxes.clear();
}

void CPUGillespieParallel::handleBoxReactions() {
    using promise_t = std::promise<std::set<data_t::index_t>>;
    using promise_new_particles_t = std::promise<data_t::update_t>;
    using promise_records = std::promise<std::vector<record_t>>;

    auto worker = [this](SlicedBox &box, ctx_t ctx, data_t *data, nl_t * nl, promise_t update,
                         promise_new_particles_t newParticles, promise_records promiseRecords) {
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
                if(ctx.recordReactionsWithPositions()) {
                    std::vector<record_t> records;
                    auto result = handleEventsGillespie(kernel, timeStep, false, approximateRate, std::move(localEvents),
                                                        &records);
                    newParticles.set_value(std::move(result));
                    promiseRecords.set_value(std::move(records));
                } else {
                    auto result = handleEventsGillespie(kernel, timeStep, false, approximateRate, std::move(localEvents),
                                                        nullptr);
                    newParticles.set_value(std::move(result));
                    std::vector<record_t> no_records;
                    promiseRecords.set_value(std::move(no_records));
                }

            }

        }
        update.set_value(std::move(problematic));
    };

    std::vector<std::future<std::set<data_t::index_t>>> updates;
    std::vector<std::future<data_t::update_t>> newParticles;
    std::vector<std::future<std::vector<record_t>>> records;
    auto &stateModel = kernel->getCPUKernelStateModel();
    {
        //readdy::util::Timer t ("\t run threads");
        std::vector<thd::scoped_thread> threads;
        for (unsigned int i = 0; i < kernel->getNThreads(); ++i) {
            // nboxes == nthreads
            promise_t promise;
            updates.push_back(promise.get_future());
            promise_new_particles_t promiseParticles;
            newParticles.push_back(promiseParticles.get_future());
            promise_records promiseRecords;
            records.push_back(promiseRecords.get_future());
            threads.push_back(
                    thd::scoped_thread(std::thread(
                            worker, std::ref(boxes[i]),
                            std::ref(kernel->getKernelContext()),
                            stateModel.getParticleData(),
                            stateModel.getNeighborList(),
                            std::move(promise),
                            std::move(promiseParticles),
                            std::move(promiseRecords)
                    ))
            );
        }
    }
    std::vector<std::vector<record_t>> attainedRecords;
    attainedRecords.reserve(records.size()+1);
    {
        //readdy::util::Timer t ("\t fix marked");
        auto &data = *stateModel.getParticleData();
        const auto& d2 = kernel->getKernelContext().getDistSquaredFun();
        auto neighbor_list = stateModel.getNeighborList();
        std::vector<event_t> evilEvents{};
        double alpha = 0;
        long n_local_problematic = 0;
        for (auto&& update : updates) {
            auto local_problematic = std::move(update.get());
            n_local_problematic += local_problematic.size();
            gatherEvents(kernel, std::move(local_problematic), neighbor_list, data, alpha, evilEvents, d2);
        }
        if(kernel->getKernelContext().recordReactionsWithPositions()) {
            std::vector<record_t> newRecords;
            auto newProblemParticles = handleEventsGillespie(kernel, timeStep, false, approximateRate, std::move(evilEvents), &newRecords);
            const auto &fixPos = kernel->getKernelContext().getFixPositionFun();
            for (auto &&future : newParticles) {
                neighbor_list->updateData(std::move(future.get()));
            }
            neighbor_list->updateData(std::move(newProblemParticles));
            std::size_t totalRecordSize = newRecords.size();
            attainedRecords.push_back(std::move(newRecords));
            for(auto &&future : records) {
                attainedRecords.push_back(std::move(future.get()));
                totalRecordSize += attainedRecords.back().size();
            }
            auto& modelRecords = kernel->getCPUKernelStateModel().reactionRecords();
            modelRecords.clear();
            modelRecords.reserve(totalRecordSize);
            for(const auto& r : attainedRecords) {
                modelRecords.insert(modelRecords.end(), r.begin(), r.end());
            }
        } else {
            auto newProblemParticles = handleEventsGillespie(kernel, timeStep, false, approximateRate, std::move(evilEvents), nullptr);
            const auto &fixPos = kernel->getKernelContext().getFixPositionFun();
            for (auto &&future : newParticles) {
                neighbor_list->updateData(std::move(future.get()));
            }
            neighbor_list->updateData(std::move(newProblemParticles));
        }

    }

}

void CPUGillespieParallel::findProblematicParticles(
        data_t::index_t index, const SlicedBox &box, ctx_t ctx,
        const data_t &data, nl_t* nl, std::set<data_t::index_t> &problematic,
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

    for (const auto neighbor : nl->find_neighbors(index)) {
        const auto& neighborEntry = data.entry_at(neighbor);
        const auto &reactions = ctx.getOrder2Reactions(me.type, neighborEntry.type);
        if (!reactions.empty()) {
            const auto distSquared = d2(neighborEntry.position(), me.position());
            for (const auto &r : reactions) {
                if (r->getRate() > 0 && distSquared < r->getEductDistanceSquared()) {
                    const bool neighborProblematic = problematic.find(neighbor) != problematic.end();
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
        for (const auto x_neighbor : nl->find_neighbors(x)) {
            const auto& x_neighbor_entry = data.entry_at(x_neighbor);
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
                const auto &reactions = ctx.getOrder2Reactions(x_entry.type, x_neighbor_entry.type);
                if (!reactions.empty()) {
                    //BOOST_LOG_TRIVIAL(debug) << "\t\t neighbor had potentially conflicting reactions";
                    const bool alreadyProblematic = problematic.find(x_neighbor) != problematic.end();
                    if (alreadyProblematic) {
                        //BOOST_LOG_TRIVIAL(debug) << "\t\t neighbor was already found, ignore";
                        continue;
                    }
                    const auto distSquared = d2(x_entry.position(), x_neighbor_entry.position());
                    for (const auto &reaction : reactions) {
                        if (reaction->getRate() > 0 && distSquared < reaction->getEductDistanceSquared()) {
                            // we have a problematic particle!
                            problematic.insert(x_neighbor);
                            bfs.push(x_neighbor);
                            break;
                        }
                    }
                }
            }
        }
    }

}

double CPUGillespieParallel::getMaxReactionRadius() const {
    return maxReactionRadius;
}

double CPUGillespieParallel::getBoxWidth() const {
    return boxWidth;
}

unsigned int CPUGillespieParallel::getLongestAxis() const {
    return longestAxis;
}

unsigned int CPUGillespieParallel::getOtherAxis1() const {
    return otherAxis1;
}

unsigned int CPUGillespieParallel::getOtherAxis2() const {
    return otherAxis2;
}

void CPUGillespieParallel::setApproximateRate(bool approximateRate) {
    CPUGillespieParallel::approximateRate = approximateRate;
}

CPUGillespieParallel::~CPUGillespieParallel() = default;
}

}
}
}
}