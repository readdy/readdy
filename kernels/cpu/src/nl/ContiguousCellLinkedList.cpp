/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          *
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
 * @file ContiguousCellLinkedList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12/24/17
 */
#include "readdy/kernel/cpu/nl/ContiguousCellLinkedList.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

ContiguousCellLinkedList::ContiguousCellLinkedList(data_type &data, const readdy::model::Context &context,
                                                   thread_pool &pool)
        : CellLinkedList(data, context, pool) {}

void ContiguousCellLinkedList::fillBins(const util::PerformanceNode &node) {
    auto t = node.timeit();

    const auto nCells = _cellIndex.size();
    _blockNParticles.resize(nCells);
    std::for_each(_blockNParticles.begin(), _blockNParticles.end(), [](auto &x) { *x = 0; });

    auto maxParticlesPerCell = getMaxCounts(node.subnode("getMaxCounts"));
    if(maxParticlesPerCell == 0) {
        _binsIndex = util::Index2D(nCells, 0_z);
    }
    _binsIndex = util::Index2D(nCells, maxParticlesPerCell);
    _bins.resize(0);
    _bins.resize(_binsIndex.size());

    auto boxSize = _context.get().boxSize();

    auto particleInBox = [boxSize](const Vec3 &pos) {
        return -.5*boxSize[0] <= pos.x && .5*boxSize[0] > pos.x
               && -.5*boxSize[1] <= pos.y && .5*boxSize[1] > pos.y
               && -.5*boxSize[2] <= pos.z && .5*boxSize[2] > pos.z;
    };

    {
        auto &blockNParticles = _blockNParticles;
        const auto &data = _data.get();
        auto cellSize = _cellSize;
        const auto &cellIndex = _cellIndex;
        auto &bins = _bins;
        const auto &binsIndex = _binsIndex;
        auto worker =
                [&data, &blockNParticles, &bins, cellIndex, cellSize, boxSize, nCells, binsIndex, maxParticlesPerCell,
                        particleInBox]
                (std::size_t tid, std::size_t begin_pidx, std::size_t end_pidx) {
            std::vector<count_type> cellCounts;
            cellCounts.resize(nCells);

            auto it = data.begin() + begin_pidx;
            auto pidx = begin_pidx;
            while (it != data.begin() + end_pidx) {
                if (!it->deactivated && particleInBox(it->pos)) {
                    const auto i = static_cast<std::size_t>(std::floor((it->pos.x + .5 * boxSize[0]) / cellSize.x));
                    const auto j = static_cast<std::size_t>(std::floor((it->pos.y + .5 * boxSize[1]) / cellSize.y));
                    const auto k = static_cast<std::size_t>(std::floor((it->pos.z + .5 * boxSize[2]) / cellSize.z));
                    const auto boxIdx = cellIndex(i, j, k);
                    auto size = (*blockNParticles.at(boxIdx))++;
                    bins.at(binsIndex(boxIdx, size)) = pidx;
                }
                ++pidx;
                ++it;
            }
        };

        auto &pool = _pool.get();
        std::vector<util::thread::joining_future<void>> futures;
        futures.reserve(pool.size());
        const auto grainSize = data.size() / pool.size();
        auto it = 0_z;
        for (std::size_t i = 0; i < pool.size()-1; ++i) {
            auto itNext = it + grainSize;
            if(it != itNext) {
                futures.emplace_back(pool.push(worker, it, itNext));
            }
            it = itNext;
        }
        if(it != data.size()) futures.emplace_back(pool.push(worker, it, data.size()));
    }

}

ContiguousCellLinkedList::count_type ContiguousCellLinkedList::getMaxCounts(const util::PerformanceNode &node) {
    auto t = node.timeit();

    const auto &boxSize = _context.get().boxSize();
    const auto &data = _data.get();
    auto &pool = _pool.get();
    const auto grainSize = data.size() / pool.size();
    const auto cellSize = _cellSize;
    const auto &cellIndex = _cellIndex;
    auto nCells = cellIndex.size();

    auto particleInBox = [boxSize](const Vec3 &pos) {
        return -.5*boxSize[0] <= pos.x && .5*boxSize[0] > pos.x
               && -.5*boxSize[1] <= pos.y && .5*boxSize[1] > pos.y
               && -.5*boxSize[2] <= pos.z && .5*boxSize[2] > pos.z;
    };

    block_n_particles_type blockNParticles;
    blockNParticles.resize(nCells);
    std::for_each(blockNParticles.begin(), blockNParticles.end(), [](auto &x) { *x = 0; });

    {
        auto worker = [&data, cellIndex, cellSize, boxSize, nCells, &blockNParticles, particleInBox]
                (std::size_t tid, std::size_t begin_pidx, std::size_t end_pidx) {
            std::vector<count_type> cellCounts;
            cellCounts.resize(nCells);

            auto it = data.begin() + begin_pidx;
            while (it != data.begin() + end_pidx) {
                const auto &entry = *it;
                if (!entry.deactivated && particleInBox(entry.pos)) {
                    const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / cellSize.x));
                    const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / cellSize.y));
                    const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / cellSize.z));
                    (*blockNParticles.at(cellIndex(i, j, k)))++;
                }
                ++it;
            }

        };

        std::vector<util::thread::joining_future<void>> futures;
        futures.reserve(pool.size());
        auto it = 0_z;
        for (std::size_t i = 0; i < pool.size()-1; ++i) {
            auto itNext = it + grainSize;
            if(it != itNext) {
                futures.emplace_back(pool.push(worker, it, itNext));
            }
            it = itNext;
        }
        if(it != data.size()) {
            futures.emplace_back(pool.push(worker, it, data.size()));
        }
    }

    auto maxCounts = *std::max_element(blockNParticles.begin(), blockNParticles.end(),
                                       [](const auto& n1, const auto &n2) -> bool { return *n1 < *n2; });
    log::debug("found cell with {} particles", *maxCounts);
    return *maxCounts;
}

}
}
}
}
