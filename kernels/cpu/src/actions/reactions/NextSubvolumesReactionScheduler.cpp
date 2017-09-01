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
 * @file SingleCPUNextSubvolumes.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.09.16
 */

#include <readdy/kernel/cpu/actions/reactions/ReactionUtils.h>
#include <readdy/common/numeric.h>
#include "readdy/kernel/cpu/actions/reactions/NextSubvolumesReactionScheduler.h"

using particle_type = readdy::model::Particle;
namespace rnd = readdy::model::rnd;

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

struct CPUNextSubvolumes::ReactionEvent {
    int order;
    unsigned int type1, type2;
    std::size_t reactionIndex;
    scalar reactionRate;
    scalar cumulativeRate {0};

    explicit ReactionEvent(const int order, const size_t reactionIndex, const scalar reactionRate, unsigned int t1,
                           unsigned int t2)
            : order(order), reactionIndex(reactionIndex), reactionRate(reactionRate), type1(t1), type2(t2) {}

    ReactionEvent() : ReactionEvent(-1, 0, 0, 0, 0) {}

    bool isValid() const {
        return order >= 0;
    }
};

struct CPUNextSubvolumes::GridCell {
    using particle_index = data_t::index_t;

    cell_index_t i, j, k, id;
    std::vector<const GridCell *> neighbors;
    std::unordered_map<particle_type::type_type, std::vector<particle_index>> particles;
    std::vector<unsigned long> typeCounts;
    scalar cellRate;
    scalar timestamp;
    ReactionEvent nextEvent{};

    GridCell(const cell_index_t i_, const cell_index_t j_, const cell_index_t k_, const cell_index_t id_,
             const unsigned long nTypes)
            : i(i_), j(j_), k(k_), id(id_), neighbors(27), particles({}), cellRate(0), typeCounts(nTypes), timestamp(0) {
    }

    void addNeighbor(GridCell const *const cell) {
        neighbors.push_back(cell);
    }
};

void CPUNextSubvolumes::perform(bool measure, const std::string &measureLabel) {
    util::RAIITimer timer(measure, measureLabel);
    eventQueue.clear();
    setUpGrid();
    assignParticles();
    setUpEventQueue();
    evaluateReactions();
}

CPUNextSubvolumes::CPUNextSubvolumes(const CPUKernel *const kernel, scalar timeStep)
        : super(timeStep), kernel(kernel), cells({}), nCells({}), cellSize(), eventQueue({}) {}

void CPUNextSubvolumes::setUpGrid() {
    if (cells.empty()) {
        const auto &simBoxSize = kernel->getKernelContext().getBoxSize();
        const auto minCellWidth = getMaxReactionRadius();
        const auto nTypes = kernel->getKernelContext().particle_types().n_types();
        for (unsigned int i = 0; i < 3; ++i) {
            nCells[i] = minCellWidth > 0 ? static_cast<unsigned int>(floor(simBoxSize[i] / minCellWidth)) : 1;
            if (nCells[i] == 0) nCells[i] = 1;
            cellSize[i] = simBoxSize[i] / nCells[i];
        }
        for (cell_index_t i = 0; i < nCells[0]; ++i) {
            for (cell_index_t j = 0; j < nCells[1]; ++j) {
                for (cell_index_t k = 0; k < nCells[2]; ++k) {
                    cells.push_back({i, j, k, k + j * nCells[2] + i * nCells[2] * nCells[1], nTypes});
                }
            }
        }
        for (auto &cell : cells) {
            setUpNeighbors(cell);
        }
    }
}

void CPUNextSubvolumes::assignParticles() {
    // todo this can be easily parallelized
    const auto data = kernel->getCPUKernelStateModel().getParticleData();
    std::for_each(cells.begin(), cells.end(), [](GridCell &cell) {
        cell.particles.clear();
        cell.typeCounts.clear();
        cell.cellRate = 0;
    });

    std::size_t idx = 0;
    for (auto& e : *data) {
        if(!e.is_deactivated()) {
            auto box = getCell(e.position());
            if (box) {
                box->particles[e.type].push_back(idx);
                ++(box->typeCounts[e.type]);
            }
        }
        ++idx;
    }

}

void CPUNextSubvolumes::evaluateReactions() {
    std::vector<GridCell *>::size_type n_cells_done = 0;
    const auto& ctx = kernel->getKernelContext();
    auto data = kernel->getCPUKernelStateModel().getParticleData();
    auto neighbor_list = kernel->getCPUKernelStateModel().getNeighborList();

    const auto comparator = [](const GridCell *c1, const GridCell *c2) {
        return c1->cellRate < c2->cellRate;
    };
    const auto get_event_queue_end = [this, &n_cells_done] {
        return eventQueue.end() - n_cells_done;
    };

    while (eventQueue.begin() != get_event_queue_end()) {
        std::pop_heap(eventQueue.begin(), get_event_queue_end(), comparator);
        auto currentCell = eventQueue.back();
        bool performedSomething = false;
        if (currentCell->timestamp < timeStep) {
            const auto& particles = currentCell->particles;
            const auto& event = currentCell->nextEvent;
            if(event.isValid()) {
                const auto findType1 = particles.find(event.type1);
                if(findType1 != particles.end()) {
                    const auto& particlesType1 = findType1->second;
                    if (!particlesType1.empty()) {
                        //const auto pidx_it = particlesType1.begin() +
                        //        static_cast<std::size_t>(rnd::uniform_int(0L, particlesType1.size()));
                        //const auto particle_idx = *pidx_it;
                        // todo execute the reaction, and update cells and event queue (see slides)
                        switch (event.order) {
                            case 1: {/*
                                // todo update neighbor list with this particle
                                scpu_data::new_entries newParticles{};
                                std::vector<scpu_data::entry_index> decayedEntries {};

                                auto reaction = ctx.getOrder1Reactions(event.type1)[event.reactionIndex];
                                performReaction(*data, particle_idx, 0, newParticles, decayedEntries, reaction);
                                performedSomething = true;
                                --currentCell->typeCounts[event.type1];
                                currentCell->particles[event.type1].erase(pidx_it);
                                if (reaction->getNProducts() > 0) {
                                    auto c1Out = getCell(p1->pos);
                                    ++c1Out->typeCounts[p1->type];
                                    c1Out->particles[p1->type].push_back(p1);
                                    if (reaction->getNProducts() == 2) {
                                        const auto& entry2 = newEntries.back();
                                        auto c2Out = getCell(entry2.pos);
                                        auto idx = data->addEntry(entry2);
                                        neighbor_list->insert(*data, idx);
                                        ++c2Out->typeCounts[entry2.type];
                                        c2Out->particles[entry2.type].push_back(idx);
                                    }
                                } else {
                                    data->removeEntry(p1);
                                    neighbor_list->remove(p1);
                                }
                                break;
                            */}
                            case 2: {
                                /*auto reaction = ctx.getOrder2Reactions(event.type1, event.type2)[event.reactionIndex];
                                const auto findType2 = particles.find(event.type2);
                                if(findType2 != particles.end()) {
                                    const auto& particlesType2 = findType2->second;
                                    const auto p2_it = rnd::random_element(particlesType2.begin(), particlesType2.end());
                                    const auto p2_idx = *p2_it;
                                    const auto p2 = (*data)[p2_idx];

                                    reaction->perform(p1, p2, pOut1, pOut2);
                                    performedSomething = true;

                                    auto c1Out = getCell(pOut1.getPos());
                                    ++c1Out->typeCounts[pOut1.getType()];
                                    c1Out->particles[pOut1.getType()].push_back(p2_idx);
                                    if (reaction->getNProducts() == 2) {
                                        auto c2Out = getCell(pOut2.getPos());
                                        data->addParticle(pOut2);
                                        neighbor_list->insert(*data, data->size()-1);
                                        ++c2Out->typeCounts[pOut2.getType()];
                                        c2Out->particles[pOut2.getType()].push_back(data->size()-1);
                                    }

                                }*/
                                break;
                            }
                            default: {
                                log::error("encountered event with order > 2! (order was {})",
                                                      event.order);
                            }
                        }
                    }
                } else {
                    log::debug("did not find any particle of type {} in current box", event.type1);
                }
            } else {
                log::error("The event was not set previously, should not happen! (?)");
            }

            if(performedSomething) {
                // todo do something
            }

        } else {
            ++n_cells_done;
        }

        std::push_heap(eventQueue.begin(), get_event_queue_end(), comparator);
    }
}

scalar CPUNextSubvolumes::getMaxReactionRadius() const {
    scalar maxReactionRadius = 0.0;
    for (auto &&e : kernel->getKernelContext().reactions().order2_flat()) {
        maxReactionRadius = std::max(maxReactionRadius, e->getEductDistance());
    }
    return maxReactionRadius;
}

void CPUNextSubvolumes::setUpNeighbors(CPUNextSubvolumes::GridCell &cell) {
    for (int _i = -1; _i < 2; ++_i) {
        for (int _j = -1; _j < 2; ++_j) {
            for (int _k = -1; _k < 2; ++_k) {
                // i am not my own neighbor
                if (!(_i == 0 && _j == 0 && _k == 0)) {
                    cell.addNeighbor(getCell(cell.i + _i, cell.j + _j, cell.k + _k));
                }
            }
        }
    }
}

CPUNextSubvolumes::GridCell *
CPUNextSubvolumes::getCell(signed_cell_index_t i, signed_cell_index_t j, signed_cell_index_t k) {
    const auto &periodic = kernel->getKernelContext().getPeriodicBoundary();
    if (periodic[0]) i = static_cast<cell_index_t>(readdy::util::numeric::positive_modulo(i, nCells[0]));
    else if (i < 0 || i >= nCells[0]) return nullptr;
    if (periodic[1]) j = static_cast<cell_index_t>(readdy::util::numeric::positive_modulo(j, nCells[1]));
    else if (j < 0 || j >= nCells[1]) return nullptr;
    if (periodic[2]) k = static_cast<cell_index_t>(readdy::util::numeric::positive_modulo(k, nCells[2]));
    else if (k < 0 || k >= nCells[2]) return nullptr;
    return &*(cells.begin() + k + j * nCells[2] + i * nCells[2] * nCells[1]);
}

void CPUNextSubvolumes::setUpEventQueue() {
    const auto &ctx = kernel->getKernelContext();
    const auto dt = timeStep;
    {
        // todo: this can easily be parallelized
        for (auto &cell : cells) {
            setUpCell(cell);
        }
    }
    {
        // fill up event queue with viable cells
        for (auto &cell : cells) {
            if (cell.timestamp < dt) eventQueue.push_back(&cell);
        }
    }

    std::for_each(cells.begin(), cells.end(), [this](GridCell &cell) { eventQueue.push_back(&cell); });
    const auto comparator = [](const GridCell *c1, const GridCell *c2) {
        return c1->timestamp > c2->timestamp;
    };
    std::make_heap(eventQueue.begin(), eventQueue.end(), comparator);
}

void CPUNextSubvolumes::setUpCell(CPUNextSubvolumes::GridCell &cell) {
    const auto &ctx = kernel->getKernelContext();
    std::vector<ReactionEvent> events;
    // order 1 reactions
    for (const auto& it : cell.particles) {
        const auto pType = it.first;
        for(auto i = 0; i < it.second.size(); ++i) {
            std::size_t reactionIdx = 0;
            for (const auto reactionOrder1 : ctx.reactions().order1_by_type(pType)) {
                const auto rateUpdate = cell.typeCounts[pType] * reactionOrder1->getRate();
                if (rateUpdate > 0) {
                    cell.cellRate += rateUpdate;
                    auto evt = ReactionEvent(1, reactionIdx, rateUpdate, pType, 0);
                    evt.cumulativeRate = evt.reactionRate;
                    if (!events.empty()) evt.cumulativeRate += events.back().cumulativeRate;
                    events.push_back(std::move(evt));
                }
                ++reactionIdx;
            }
        }
    }
    // order 2 reactions
    {
        std::size_t reactionIdx = 0;
        for (const auto reactionOrder2 : ctx.reactions().order2_flat()) {
            for (int perm = 0; perm < 2; ++perm) {
                const auto typeA = reactionOrder2->getEducts()[(0 + perm) % 2];
                const auto typeB = reactionOrder2->getEducts()[(1 + perm) % 2];
                const auto neighborSum = std::accumulate(
                        cell.neighbors.begin(), cell.neighbors.end(), 0,
                        [typeB](unsigned long acc, const GridCell *const neighborCell) {
                            return acc + neighborCell->typeCounts[typeB];
                        }
                );
                const auto rateUpdate = .5 * static_cast<scalar>(cell.typeCounts[typeA] * neighborSum);
                if (rateUpdate > 0) {
                    auto evt = ReactionEvent(1, reactionIdx, rateUpdate, typeA, typeB);
                    cell.cellRate += rateUpdate;
                    evt.cumulativeRate = evt.reactionRate;
                    if (!events.empty()) evt.cumulativeRate += events.back().cumulativeRate;
                    events.push_back(std::move(evt));
                }
            }
            ++reactionIdx;
        }
    }
    cell.timestamp = rnd::exponential(cell.cellRate);
    // select next event
    {
        const auto x = rnd::uniform_real(static_cast<scalar>(0.), events.back().cumulativeRate);
        const auto eventIt = std::lower_bound(
                events.begin(), events.end(), x, [](const ReactionEvent &elem1, scalar elem2) {
                    return elem1.cumulativeRate < elem2;
                }
        );
        if (eventIt != events.end()) {
            cell.nextEvent = std::move(events[eventIt - events.begin()]);
        } else {
            log::error("next subvolumes: the next event was events.end()!");
        }
    }
}

CPUNextSubvolumes::GridCell *CPUNextSubvolumes::getCell(const readdy::model::Vec3 &particlePosition) {
    const auto& simBoxSize = kernel->getKernelContext().getBoxSize();
    const auto i = static_cast<cell_index_t>(floor((particlePosition[0] + .5*simBoxSize[0]) / cellSize[0]));
    const auto j = static_cast<cell_index_t>(floor((particlePosition[1] + .5*simBoxSize[1]) / cellSize[1]));
    const auto k = static_cast<cell_index_t>(floor((particlePosition[2] + .5*simBoxSize[2]) / cellSize[2]));
    return getCell(i, j, k);
}

CPUNextSubvolumes::~CPUNextSubvolumes() = default;


}
}
}
}
}