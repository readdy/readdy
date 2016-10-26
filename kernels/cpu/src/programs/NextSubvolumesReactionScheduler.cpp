/**
 * << detailed description >>
 *
 * @file SingleCPUNextSubvolumes.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.09.16
 */

#include "readdy/kernel/cpu/programs/NextSubvolumesReactionScheduler.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

struct NextSubvolumes::ReactionEvent {
    int order;
    unsigned int type1, type2;
    std::size_t reactionIndex;
    double reactionRate;
    double cumulativeRate;

    explicit ReactionEvent(const int order, const size_t reactionIndex, const double reactionRate, unsigned int t1,
                           unsigned int t2)
            : order(order), reactionIndex(reactionIndex), reactionRate(reactionRate), type1(t1), type2(t2) {}

    ReactionEvent() : ReactionEvent(-1, 0, 0, 0, 0) {}

    bool isValid() {
        return order >= 0;
    }
};

struct NextSubvolumes::GridCell {
    const cell_index_t i, j, k, id;
    std::vector<const GridCell *> neighbors;
    std::vector<unsigned long> particles;
    std::vector<unsigned long> typeCounts;
    double cellRate;
    double timestamp;
    ReactionEvent nextEvent {};

    GridCell(const cell_index_t i, const cell_index_t j, const cell_index_t k, const cell_index_t id,
                           const unsigned long nTypes)
            : i(i), j(j), k(k), id(id), neighbors({27}), particles({}), cellRate(0), typeCounts(nTypes), timestamp(0) {
    }

    void addNeighbor(GridCell const *const cell) {
        neighbors.push_back(cell);
    }
};

void NextSubvolumes::execute() {
    eventQueue.clear();
    setUpGrid();
    assignParticles();
    setUpEventQueue();
    evaluateReactions();
}

NextSubvolumes::NextSubvolumes(const CPUKernel *const kernel)
        : kernel(kernel), cells({}), nCells({}), cellSize(), eventQueue({}) {}

void NextSubvolumes::setUpGrid() {
    if (cells.empty()) {
        const auto &simBoxSize = kernel->getKernelContext().getBoxSize();
        const auto minCellWidth = getMaxReactionRadius();
        const auto nTypes = kernel->getKernelContext().getAllRegisteredParticleTypes().size();
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

void NextSubvolumes::assignParticles() {
    const auto simBoxSize = kernel->getKernelContext().getBoxSize();
    const auto data = kernel->getKernelStateModel().getParticleData();
    std::for_each(cells.begin(), cells.end(), [](GridCell &cell) {
        cell.particles.clear();
        cell.typeCounts.clear();
        cell.cellRate = 0;
    });

    const auto shift = .5 * readdy::model::Vec3(simBoxSize[0], simBoxSize[1], simBoxSize[2]);
    unsigned long idx = 0;
    auto it_type = data->begin_types();
    for (auto it = data->begin_positions(); it != data->end_positions(); ++it) {
        const auto pos_shifted = *it + shift;
        const auto i = static_cast<cell_index_t>(floor(pos_shifted[0] / cellSize[0]));
        const auto j = static_cast<cell_index_t>(floor(pos_shifted[1] / cellSize[1]));
        const auto k = static_cast<cell_index_t>(floor(pos_shifted[2] / cellSize[2]));
        auto box = getCell(i, j, k);
        if (box) {
            box->particles.push_back(idx);
            ++(box->typeCounts[*it_type]);
        }
        ++idx;
        ++it_type;
    }

}

void NextSubvolumes::evaluateReactions() {
    const auto comparator = [](const GridCell* c1, const GridCell* c2) {
        return c1->cellRate < c2->cellRate;
    };
    while(!eventQueue.empty()) {
        std::pop_heap(eventQueue.begin(), eventQueue.end(), comparator);
        auto currentCell = eventQueue.back();
        eventQueue.pop_back();
        {
            // todo execute the reaction, and update cells and event queue (see slides)
        }
    }
}

double NextSubvolumes::getMaxReactionRadius() const {
    double maxReactionRadius = 0.0;
    for (auto &&e : kernel->getKernelContext().getAllOrder2Reactions()) {
        maxReactionRadius = std::max(maxReactionRadius, e->getEductDistance());
    }
    return maxReactionRadius;
}

void NextSubvolumes::setUpNeighbors(NextSubvolumes::GridCell &cell) {
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

NextSubvolumes::GridCell *
NextSubvolumes::getCell(cell_index_t i, cell_index_t j, cell_index_t k) {
    const auto &periodic = kernel->getKernelContext().getPeriodicBoundary();
    if (periodic[0]) i = static_cast<cell_index_t>(readdy::util::numeric::positive_modulo(i, nCells[0]));
    else if (i < 0 || i >= nCells[0]) return nullptr;
    if (periodic[1]) j = static_cast<cell_index_t>(readdy::util::numeric::positive_modulo(j, nCells[1]));
    else if (j < 0 || j >= nCells[1]) return nullptr;
    if (periodic[2]) k = static_cast<cell_index_t>(readdy::util::numeric::positive_modulo(k, nCells[2]));
    else if (k < 0 || k >= nCells[2]) return nullptr;
    return &cells[k + j * nCells[2] + i * nCells[2] * nCells[1]];
}

void NextSubvolumes::setUpEventQueue() {
    const auto &ctx = kernel->getKernelContext();
    const auto dt = ctx.getTimeStep();
    {
        // todo: this can easily be parallelized
        for (auto &cell : cells) {
           setUpCell(cell);
        }
    }
    {
        // fill up event queue with viable cells
        for(auto &cell : cells) {
            if(cell.timestamp < dt) eventQueue.push_back(&cell);
        }
    }

    std::for_each(cells.begin(), cells.end(), [this](GridCell &cell) { eventQueue.push_back(&cell); });
    const auto comparator = [](const GridCell* c1, const GridCell* c2) {
        return c1->timestamp > c2->timestamp;
    };
    std::make_heap(eventQueue.begin(), eventQueue.end(), comparator);
}

void NextSubvolumes::setUpCell(NextSubvolumes::GridCell &cell) {
    const auto &ctx = kernel->getKernelContext();
    const auto data = kernel->getKernelStateModel().getParticleData();
    std::vector<ReactionEvent> events;
    // order 1 reactions
    for (const auto pIdx : cell.particles) {
        const auto pType = *(data->begin_types() + pIdx);
        std::size_t reactionIdx = 0;
        for (const auto reactionOrder1 : ctx.getOrder1Reactions(pType)) {
            const auto rateUpdate = cell.typeCounts[pType] * reactionOrder1->getRate();
            if(rateUpdate > 0) {
                cell.cellRate += rateUpdate;
                auto evt = ReactionEvent(1, reactionIdx, rateUpdate, pType, 0);
                evt.cumulativeRate = evt.reactionRate;
                if(!events.empty()) evt.cumulativeRate += events.back().cumulativeRate;
                events.push_back(std::move(evt));
            }
            ++reactionIdx;
        }
    }
    // order 2 reactions
    {
        std::size_t reactionIdx = 0;
        for (const auto reactionOrder2 : ctx.getAllOrder2Reactions()) {
            for (int perm = 0; perm < 2; ++perm) {
                const auto typeA = reactionOrder2->getEducts()[(0 + perm) % 2];
                const auto typeB = reactionOrder2->getEducts()[(1 + perm) % 2];
                const auto neighborSum = std::accumulate(
                        cell.neighbors.begin(), cell.neighbors.end(), 0,
                        [typeB](unsigned long acc, const GridCell *const neighborCell) {
                            return acc + neighborCell->typeCounts[typeB];
                        }
                );
                const auto rateUpdate = .5 * static_cast<double>(cell.typeCounts[typeA] * neighborSum);
                if (rateUpdate > 0) {
                    auto evt = ReactionEvent(1, reactionIdx, rateUpdate, typeA, typeB);
                    cell.cellRate += rateUpdate;
                    evt.cumulativeRate = evt.reactionRate;
                    if(!events.empty()) evt.cumulativeRate += events.back().cumulativeRate;
                    events.push_back(std::move(evt));
                }
            }
            ++reactionIdx;
        }
    }
    cell.timestamp = readdy::model::rnd::exponential(cell.cellRate);
    // select next event
    {
        const auto x = readdy::model::rnd::uniform(0, events.back().cumulativeRate);
        const auto eventIt = std::lower_bound(
                events.begin(), events.end(), x, [](const ReactionEvent &elem1, double elem2) {
                    return elem1.cumulativeRate < elem2;
                }
        );
        if(eventIt != events.end()) {
            cell.nextEvent = std::move(events[eventIt - events.begin()]);
        } else {
            log::console()->error("next subvolumes: the next event was events.end()!");
        }
    }
}

NextSubvolumes::~NextSubvolumes() = default;


}
}
}
}
}