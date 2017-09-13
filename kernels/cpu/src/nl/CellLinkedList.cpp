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

CellLinkedList::CellLinkedList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                               const readdy::util::thread::Config &config)
        : _data(data), _context(context), _config(config) {}


void CellLinkedList::setUp(scalar skin, std::uint8_t radius, const util::PerformanceNode &node) {
    if(!_is_set_up || _skin != skin || _radius != radius) {
        auto t = node.timeit();

        _skin = skin;
        _radius = radius;
        _max_cutoff = _context.get().calculateMaxCutoff();
        _max_cutoff_skin_squared = (_max_cutoff + _skin) * (_max_cutoff + _skin);
        if(_max_cutoff > 0) {
            auto size = _context.get().boxSize();
            auto desiredWidth = static_cast<scalar>((_max_cutoff + _skin)/static_cast<scalar>(radius));
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
                auto pbc = _context.get().periodicBoundaryConditions();
                {
                    int r = _radius;
                    // local adjacency
                    std::vector<std::size_t> adj;
                    adj.reserve(1 + nAdjacentCells);

                    for (int i = 0; i < _cellIndex[0]; ++i) {
                        for (int j = 0; j < _cellIndex[1]; ++j) {
                            for (int k = 0; k < _cellIndex[2]; ++k) {
                                auto cellIdx = _cellIndex(static_cast<std::size_t>(i), static_cast<std::size_t>(j),
                                                          static_cast<std::size_t>(k));

                                adj.clear();

                                for (int ii = i - r; ii <= i + r; ++ii) {
                                    for (int jj = j - r; jj <= j + r; ++jj) {
                                        for (int kk = k - r; kk <= k + r; ++kk) {
                                            if (ii == i && jj == j && kk == k) continue;
                                            auto adj_x = ii;
                                            if (pbc[0] && adj_x < 0)
                                                adj_x = static_cast<int>(util::numeric::positive_modulo(adj_x,
                                                                                                        _cellIndex[0]));
                                            auto adj_y = jj;
                                            if (pbc[1] && adj_y < 0)
                                                adj_y = static_cast<int>(util::numeric::positive_modulo(adj_y,
                                                                                                        _cellIndex[1]));
                                            auto adj_z = kk;
                                            if (pbc[2] && adj_z < 0)
                                                adj_z = static_cast<int>(util::numeric::positive_modulo(adj_z,
                                                                                                        _cellIndex[2]));

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

const util::Index3D &CellLinkedList::cellIndex() const {
    return _cellIndex;
}

const util::Index2D &CellLinkedList::neighborIndex() const {
    return _cellNeighbors;
}

const std::vector<std::size_t> &CellLinkedList::neighbors() const {
    return _cellNeighborsContent;
}

std::size_t *CellLinkedList::neighborsBegin(std::size_t cellIndex) {
    auto beginIdx = _cellNeighbors(cellIndex, 1_z);
    return &_cellNeighborsContent.at(beginIdx);
}

const std::size_t *CellLinkedList::neighborsBegin(std::size_t cellIndex) const {
    auto beginIdx = _cellNeighbors(cellIndex, 1_z);
    return &_cellNeighborsContent.at(beginIdx);
}

std::size_t *CellLinkedList::neighborsEnd(std::size_t cellIndex) {
    return neighborsBegin(cellIndex) + nNeighbors(cellIndex);
}

const std::size_t *CellLinkedList::neighborsEnd(std::size_t cellIndex) const {
    return neighborsBegin(cellIndex) + nNeighbors(cellIndex);
}

std::size_t CellLinkedList::nNeighbors(std::size_t cellIndex) const {
    return _cellNeighborsContent.at(_cellNeighbors(cellIndex, 0_z));
}


ContiguousCellLinkedList::ContiguousCellLinkedList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                                                   const util::thread::Config &config)
        : CellLinkedList(data, context, config) {}

void ContiguousCellLinkedList::initializeBinsStructure() {
    _binsIndex = util::Index2D(_cellIndex.size(), 1+_maxParticlesPerCell);
    _bins.clear();
    _bins.resize(_binsIndex.size());
}

void ContiguousCellLinkedList::setUpBins(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if(_maxParticlesPerCell == 0) {
        // assume uniform distribution
        auto nParticles = static_cast<scalar>(_data.get().size() - _data.get().getNDeactivated());
        _maxParticlesPerCell = std::max(1_z, static_cast<std::size_t>(std::ceil(1.5 * nParticles / _context.get().boxVolume())));
    }
    initializeBinsStructure();
    fillBins(node.subnode("fillBins"));
}

void ContiguousCellLinkedList::fillBins(const util::PerformanceNode &node) {
    auto t = node.timeit();
    auto boxSize = _context.get().boxSize();
    std::size_t maxParticlesPerCell = 0;
    do {
        if(maxParticlesPerCell != 0)  {
            _maxParticlesPerCell = maxParticlesPerCell;
            initializeBinsStructure();
        }
        bool abort = false;
        std::size_t pidx = 0;
        for (const auto &entry : _data.get()) {
            if (!entry.deactivated) {
                const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / _cellSize.x));
                const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / _cellSize.y));
                const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / _cellSize.z));
                const auto boxIdx = _cellIndex(i, j, k);
                auto &size = _bins.at(_binsIndex(boxIdx, 0_z));
                size++;
                if (size > maxParticlesPerCell) {
                    maxParticlesPerCell = size;
                    abort = maxParticlesPerCell > _maxParticlesPerCell;
                }
                if (!abort) {
                    _bins.at(_binsIndex(boxIdx, size)) = pidx;
                }
            }
            ++pidx;
        }
        if(abort) {
            log::info("increasing max particles per cell to {}", maxParticlesPerCell);
        }
    } while(maxParticlesPerCell > _maxParticlesPerCell);
}

void ContiguousCellLinkedList::update(const util::PerformanceNode &node) {
    if(_max_cutoff > 0) {
        auto t = node.timeit();
        for (std::size_t box = 0; box < _cellIndex.size(); ++box) {
            _bins.at(_binsIndex(box, 0_z)) = 0;
        }
        fillBins(node.subnode("fillBins"));
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
    const auto beginIndex = _binsIndex(cellIndex, 1_z);
    return &_bins.at(beginIndex);
}

const std::size_t *ContiguousCellLinkedList::particlesBegin(std::size_t cellIndex) const {
    const auto beginIndex = _binsIndex(cellIndex, 1_z);
    return &_bins.at(beginIndex);
}

std::size_t *ContiguousCellLinkedList::particlesEnd(std::size_t cellIndex) {
    return particlesBegin(cellIndex) + nParticles(cellIndex);
}

const std::size_t *ContiguousCellLinkedList::particlesEnd(std::size_t cellIndex) const {
    return particlesBegin(cellIndex) + nParticles(cellIndex);
}

std::size_t ContiguousCellLinkedList::nParticles(std::size_t cellIndex) const {
    return _bins.at(_binsIndex(cellIndex, 0_z));
}


DynamicCellLinkedList::DynamicCellLinkedList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                                             const util::thread::Config &config)
        : CellLinkedList(data, context, config) {}

void DynamicCellLinkedList::setUpBins(const util::PerformanceNode &node) {
    auto t = node.timeit();
    {
        auto tt = node.subnode("allocate").timeit();
        _bins.clear();
        _bins.resize(_cellIndex.size());
        auto nParticles = static_cast<scalar>(_data.get().size() - _data.get().getNDeactivated());
        auto estimatedParticlesPerCell = static_cast<std::size_t>(std::ceil(
                1.5 * nParticles / _context.get().boxVolume()));
        std::for_each(_bins.begin(), _bins.end(), [estimatedParticlesPerCell](std::vector<std::size_t> &bin) {
            bin.reserve(estimatedParticlesPerCell);
        });
    }
    fillBins(node.subnode("fillBins"));
}

void DynamicCellLinkedList::fillBins(const util::PerformanceNode &node) {
    auto t = node.timeit();
    auto boxSize = _context.get().boxSize();
    std::size_t pidx = 0;
    for (const auto &entry : _data.get()) {
        if (!entry.deactivated) {
            const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / _cellSize.x));
            const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / _cellSize.y));
            const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / _cellSize.z));
            const auto boxIdx = _cellIndex(i, j, k);
            _bins.at(boxIdx).push_back(pidx);
        }
        ++pidx;
    }
}

void DynamicCellLinkedList::update(const util::PerformanceNode &node) {
    if(_max_cutoff > 0) {
        auto t = node.timeit();
        {
            auto tt = node.subnode("clear").timeit();
            std::for_each(_bins.begin(), _bins.end(), [](std::vector<std::size_t> &bin) {
                bin.clear();
            });
        }
        fillBins(node.subnode("fillBins"));
    }
}

void DynamicCellLinkedList::clear() {
    _bins.clear();
}

std::size_t *DynamicCellLinkedList::particlesBegin(std::size_t cellIndex) {
    return &*_bins.at(cellIndex).begin();
}

const std::size_t *DynamicCellLinkedList::particlesBegin(std::size_t cellIndex) const {
    return &*_bins.at(cellIndex).begin();
}

std::size_t *DynamicCellLinkedList::particlesEnd(std::size_t cellIndex) {
    return &*_bins.at(cellIndex).end();
}

const std::size_t *DynamicCellLinkedList::particlesEnd(std::size_t cellIndex) const {
    return &*_bins.at(cellIndex).end();
}

std::size_t DynamicCellLinkedList::nParticles(std::size_t cellIndex) const {
    return _bins.at(cellIndex).size();
}


}
}
}
}

