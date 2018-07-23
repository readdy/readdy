/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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
