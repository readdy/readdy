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
                                                   const util::thread::Config &config)
        : CellLinkedList(data, context, config) {}

void ContiguousCellLinkedList::fillBins(const util::PerformanceNode &node) {
    auto t = node.timeit();

    const auto nCells = _cellIndex.size();
    _blockNParticles.resize(nCells);

    auto maxParticlesPerCell = getMaxCounts(node.subnode("getMaxCounts"));
    _binsIndex = util::Index2D(nCells, maxParticlesPerCell);
    _bins.resize(0);
    _bins.resize(_binsIndex.size());


    auto boxSize = _context.get().boxSize();
    std::size_t pidx = 0;
    for (const auto &entry : _data.get()) {
        if (!entry.deactivated) {
            const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / _cellSize.x));
            const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / _cellSize.y));
            const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / _cellSize.z));
            const auto boxIdx = _cellIndex(i, j, k);
            auto size = (*_blockNParticles.at(boxIdx)).load();
            _bins.at(_binsIndex(boxIdx, size)) = pidx;
            ++(*_blockNParticles.at(boxIdx));
        }
        ++pidx;
    }
}

void ContiguousCellLinkedList::update(const util::PerformanceNode &node) {
    if(_max_cutoff > 0) {
        fillBins(node);
    }
}

void ContiguousCellLinkedList::clear() {
    _bins.clear();
}

const std::vector<std::size_t> &ContiguousCellLinkedList::bins() const {
    return _bins;
}

const util::Index2D &ContiguousCellLinkedList::binsIndex() const {
    return _binsIndex;
}

std::size_t *ContiguousCellLinkedList::particlesBegin(std::size_t cellIndex) {
    const auto beginIndex = _binsIndex(cellIndex, 0_z);
    return &_bins.at(beginIndex);
}

const std::size_t *ContiguousCellLinkedList::particlesBegin(std::size_t cellIndex) const {
    const auto beginIndex = _binsIndex(cellIndex, 0_z);
    return &_bins.at(beginIndex);
}

std::size_t *ContiguousCellLinkedList::particlesEnd(std::size_t cellIndex) {
    return particlesBegin(cellIndex) + nParticles(cellIndex);
}

const std::size_t *ContiguousCellLinkedList::particlesEnd(std::size_t cellIndex) const {
    return particlesBegin(cellIndex) + nParticles(cellIndex);
}

std::size_t ContiguousCellLinkedList::nParticles(std::size_t cellIndex) const {
    return (*_blockNParticles.at(cellIndex)).load();
}

ContiguousCellLinkedList::count_type ContiguousCellLinkedList::getMaxCounts(const util::PerformanceNode &node) {
    auto t = node.timeit();

    const auto &boxSize = _context.get().boxSize();
    const auto &data = _data.get();
    const auto grainSize = data.size() / _config.get().nThreads();
    const auto &executor = *_config.get().executor();
    const auto cellSize = _cellSize;
    const auto &cellIndex = _cellIndex;
    auto nCells = cellIndex.size();


    std::vector<count_type> blockCounts;
    {
        // initialize blockCounts
        blockCounts.resize(_config.get().nThreads());
    }

    {
        auto worker = [&data, &cellIndex, cellSize, boxSize, nCells]
                (std::size_t tid, std::size_t begin_pidx, std::size_t end_pidx, count_type &maxCount) {
            std::vector<count_type> cellCounts;
            cellCounts.resize(nCells);

            count_type localMax {0};

            auto it = data.begin() + begin_pidx - 1;
            auto pidx = begin_pidx;
            while (it != data.begin() + end_pidx - 1) {
                const auto &entry = *it;
                if (!entry.deactivated) {
                    const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / cellSize.x));
                    const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / cellSize.y));
                    const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / cellSize.z));
                    const auto cix = cellIndex(i, j, k);
                    ++cellCounts.at(cix);
                    localMax = std::max(localMax, cellCounts.at(cix));
                }
                ++pidx;
                ++it;
            }

            maxCount = localMax;
        };

        std::vector<std::function<void(std::size_t)>> executables;
        executables.reserve(_config.get().nThreads());
        auto it = 1_z;
        for (std::size_t i = 0; i < _config.get().nThreads() - 1; ++i) {
            executables.push_back(executor.pack(worker, it, it + grainSize, std::ref(blockCounts.at(i))));
            it += grainSize;
        }
        executables.push_back(executor.pack(worker, it, data.size() + 1, std::ref(blockCounts.back())));
        executor.execute_and_wait(std::move(executables));
    }

    auto maxCounts = *std::max_element(blockCounts.begin(), blockCounts.end());
    log::debug("found cell with {} particles", maxCounts);
    return maxCounts;
}

}
}
}
}
