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
 * @file CellLinkedList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/cpu/nl/CellLinkedList.h>
#include <readdy/common/numeric.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

CellLinkedList::CellLinkedList(data_type &data, const readdy::model::Context &context,
                               thread_pool &pool)
        : _data(data), _context(context), _pool(pool) {}


void CellLinkedList::setUp(scalar skin, cell_radius_type radius, const util::PerformanceNode &node) {
    if (!_is_set_up || _skin != skin || _radius != radius) {
        auto t = node.timeit();

        _skin = skin;
        _radius = radius;
        _max_cutoff = _context.get().calculateMaxCutoff();
        if (_max_cutoff > 0) {
            auto size = _context.get().boxSize();
            auto desiredWidth = static_cast<scalar>((_max_cutoff + _skin) / static_cast<scalar>(radius));
            std::array<std::size_t, 3> dims{};
            for (int i = 0; i < 3; ++i) {
                dims[i] = static_cast<unsigned int>(std::max(1., std::floor(size[i] / desiredWidth)));
                _cellSize[i] = size[i] / static_cast<scalar>(dims[i]);
            }

            _cellIndex = util::Index3D(dims[0], dims[1], dims[2]);

            {
                // set up cell adjacency list
                auto t2 = node.subnode("setUpCellNeighbors").timeit();
                std::array<std::size_t, 3> nNeighbors{{_cellIndex[0], _cellIndex[1], _cellIndex[2]}};
                for (int i = 0; i < 3; ++i) {
                    nNeighbors[i] = std::min(nNeighbors[i], static_cast<std::size_t>(2 * radius + 1));
                }
                auto nAdjacentCells = nNeighbors[0] * nNeighbors[1] * nNeighbors[2];
                _cellNeighbors = util::Index2D(_cellIndex.size(), 1 + nAdjacentCells);
                _cellNeighborsContent.resize(_cellNeighbors.size());
                {
                    auto pbc = _context.get().periodicBoundaryConditions();
                    int r = _radius;
                    // local adjacency
                    std::vector<std::size_t> adj;
                    adj.reserve(1 + nAdjacentCells);
                    for (int i = 0; i < _cellIndex[0]; ++i) {
                        for (int j = 0; j < _cellIndex[1]; ++j) {
                            for (int k = 0; k < _cellIndex[2]; ++k) {
                                auto cellIdx = _cellIndex(static_cast<std::size_t>(i), static_cast<std::size_t>(j),
                                                          static_cast<std::size_t>(k));

                                adj.resize(0);

                                for (int ii = i - r; ii <= i + r; ++ii) {
                                    for (int jj = j - r; jj <= j + r; ++jj) {
                                        for (int kk = k - r; kk <= k + r; ++kk) {
                                            if (ii == i && jj == j && kk == k) continue;
                                            auto adj_x = ii;
                                            if (pbc[0]) {
                                                auto cix = static_cast<int>(_cellIndex[0]);
                                                adj_x = (adj_x % cix + cix) % cix;
                                            }
                                            auto adj_y = jj;
                                            if (pbc[1]) {
                                                auto cix = static_cast<int>(_cellIndex[1]);
                                                adj_y = (adj_y % cix + cix) % cix;
                                            }
                                            auto adj_z = kk;
                                            if (pbc[2]) {
                                                auto cix = static_cast<int>(_cellIndex[2]);
                                                adj_z = (adj_z % cix + cix) % cix;
                                            }

                                            if (adj_x >= 0 && adj_y >= 0 && adj_z >= 0
                                                && adj_x < _cellIndex[0] && adj_y < _cellIndex[1] &&
                                                adj_z < _cellIndex[2]) {
                                                adj.push_back(_cellIndex(adj_x, adj_y, adj_z));
                                            }
                                        }
                                    }
                                }

                                std::sort(adj.begin(), adj.end());
                                adj.erase(std::unique(std::begin(adj), std::end(adj)), std::end(adj));
                                adj.erase(std::remove(std::begin(adj), std::end(adj), _cellIndex(i, j, k)),
                                          std::end(adj));

                                auto begin = _cellNeighbors(cellIdx, 0_z);
                                _cellNeighborsContent[begin] = adj.size();
                                std::copy(adj.begin(), adj.end(), &_cellNeighborsContent.at(begin + 1));
                            }
                        }
                    }
                }
            }

            if (_max_cutoff > 0) {
                setUpBins(node.subnode("setUpBins"));
            }

            _is_set_up = true;
        }
    }
}

CompactCellLinkedList::CompactCellLinkedList(data_type &data, const readdy::model::Context &context,
                                             thread_pool &pool) : CellLinkedList(data, context, pool) {}

template<>
void CompactCellLinkedList::fillBins<true>(const util::PerformanceNode &node) {
    auto t = node.timeit();
    const auto &boxSize = _context.get().boxSize();

    auto particleInBox = [boxSize](const Vec3 &pos) {
        return -.5*boxSize[0] <= pos.x && .5*boxSize[0] > pos.x
               && -.5*boxSize[1] <= pos.y && .5*boxSize[1] > pos.y
               && -.5*boxSize[2] <= pos.z && .5*boxSize[2] > pos.z;
    };

    std::size_t pidx = 1;
    for (const auto &entry : _data.get()) {
        if (!entry.deactivated && particleInBox(entry.pos)) {
            const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / _cellSize.x));
            const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / _cellSize.y));
            const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / _cellSize.z));
            const auto cellIndex = _cellIndex(i, j, k);
            _list[pidx] = *_head.at(cellIndex);
            *_head[cellIndex] = pidx;
        }
        ++pidx;
    }
}

template<>
void CompactCellLinkedList::fillBins<false>(const util::PerformanceNode &node) {
    auto t = node.timeit();

    {/*
        auto thilb = node.subnode("hilbert").timeit();
        static int i = 0;
        if(i % 100 == 0) {
            _data.get().hilbertSort(_max_cutoff + _skin);
            i = 0;
        }
        ++i;
    */}

    const auto &boxSize = _context.get().boxSize();
    const auto &data = _data.get();
    const auto grainSize = data.size() / _pool.get().size();
    const auto cellSize = _cellSize;
    const auto &cellIndex = _cellIndex;

    auto particleInBox = [boxSize](const Vec3 &pos) {
        return -.5*boxSize[0] <= pos.x && .5*boxSize[0] > pos.x
               && -.5*boxSize[1] <= pos.y && .5*boxSize[1] > pos.y
               && -.5*boxSize[2] <= pos.z && .5*boxSize[2] > pos.z;
    };

    auto &list = _list;
    auto &head = _head;

    auto worker = [&data, &cellIndex, cellSize, &list, &head, boxSize, particleInBox]
            (std::size_t tid, std::size_t begin_pidx, std::size_t end_pidx) {
        auto it = data.begin() + begin_pidx - 1;
        auto pidx = begin_pidx;
        while (it != data.begin() + end_pidx - 1) {
            const auto &entry = *it;
            if (!entry.deactivated && particleInBox(entry.pos)) {
                const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / cellSize.x));
                const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / cellSize.y));
                const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / cellSize.z));
                const auto cix = cellIndex(i, j, k);
                auto &atomic = *head.at(cix);
                // perform CAS
                auto currentHead = atomic.load();
                while (!atomic.compare_exchange_weak(currentHead, pidx)) {}
                list[pidx] = currentHead;
            }
            ++pidx;
            ++it;
        }
    };

    std::vector<util::thread::joining_future<void>> futures;
    futures.reserve(_pool.get().size());
    auto it = 1_z;
    for (auto i = 0_z; i < _pool.get().size() - 1; ++i) {
        auto itNext = it + grainSize;
        if(it != itNext) futures.emplace_back(_pool.get().push(worker, it, itNext));
        it = itNext;
    }
    futures.emplace_back(_pool.get().push(worker, it, _data.get().size()+1));
}

void CompactCellLinkedList::setUpBins(const util::PerformanceNode &node) {
    if (_max_cutoff > 0) {
        auto t = node.timeit();
        {
            auto tt = node.subnode("allocate").timeit();
            auto nParticles = _data.get().size();
            _head.clear();
            _head.resize(_cellIndex.size());
            _list.resize(0);
            _list.resize(nParticles + 1);
        }
        if (_serial) {
            fillBins<true>(node.subnode("fillBins serial"));
        } else {
            fillBins<false>(node.subnode("fillBins parallel"));
        }
    }
}

}
}
}
}
