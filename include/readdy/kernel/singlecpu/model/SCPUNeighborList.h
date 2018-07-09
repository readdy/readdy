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
 * @file SingleCPUNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#pragma once

#include <unordered_set>
#include <readdy/model/Context.h>
#include "SCPUParticleData.h"
#include <readdy/common/numeric.h>
#include <readdy/common/Index.h>
#include <readdy/common/Timer.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace model {

class BoxIterator;

class CellLinkedList {
public:
    using data_type = readdy::kernel::scpu::model::SCPUParticleData;
    using cell_radius_type = std::uint8_t;
    using HEAD = std::vector<std::size_t>;
    using LIST = std::vector<std::size_t>;
    using entry_cref = const data_type::entry_type &;
    using pair_callback = std::function<void(entry_cref, entry_cref)>;

    using iterator_bounds = std::tuple<std::size_t, std::size_t>;

    CellLinkedList(data_type &data, const readdy::model::Context &context)
            : _data(data), _context(context), _head{}, _list{}, _radius{0} {};

    void setUp(scalar skin, cell_radius_type radius, const util::PerformanceNode &node) {
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
                        auto fixBoxIx = [&](auto cix, auto boxIx, std::uint8_t axis) {
                            auto nCells = static_cast<int>(_cellIndex[axis]);
                            if (pbc[axis] && nCells > 2) {
                                return (boxIx % nCells + nCells) % nCells;
                            }
                            return boxIx;
                        };

                        int r = _radius;
                        // local adjacency
                        std::vector<std::size_t> adj;
                        adj.reserve(1 + nAdjacentCells);
                        for (int i = 0; i < _cellIndex[0]; ++i) {
                            for (int j = 0; j < _cellIndex[1]; ++j) {
                                for (int k = 0; k < _cellIndex[2]; ++k) {
                                    auto cellIdx = _cellIndex(static_cast<std::size_t>(i), static_cast<std::size_t>(j),
                                                              static_cast<std::size_t>(k));
                                    std::array<std::size_t, 3> ijk {{static_cast<std::size_t>(i), static_cast<std::size_t>(j),static_cast<std::size_t>(k)}};

                                    adj.resize(0);
                                    {
                                        std::vector<std::array<int, 3>> boxCoords{
                                                {{i + 0, j + 0, k + 1}},

                                                {{i + 0, j + 1, k - 1}},
                                                {{i + 0, j + 1, k - 0}},
                                                {{i + 0, j + 1, k + 1}},

                                                {{i + 1, j - 1, k - 1}},
                                                {{i + 1, j - 1, k + 0}},
                                                {{i + 1, j - 1, k + 1}},

                                                {{i + 1, j + 0, k - 1}},
                                                {{i + 1, j + 0, k + 0}},
                                                {{i + 1, j + 0, k + 1}},

                                                {{i + 1, j + 1, k - 1}},
                                                {{i + 1, j + 1, k + 0}},
                                                {{i + 1, j + 1, k + 1}},
                                        };

                                        std::transform(boxCoords.begin(), boxCoords.end(), boxCoords.begin(),
                                                       [&](auto arr) {
                                                           for (std::uint8_t d = 0; d < 3; ++d) {
                                                               arr.at(d) = fixBoxIx(ijk[d], arr.at(d), d);
                                                           }
                                                           return arr;
                                                       }
                                        );

                                        for (auto boxCoord : boxCoords) {
                                            if (boxCoord[0] >= 0 && boxCoord[1] >= 0 && boxCoord[2] >= 0
                                                && boxCoord[0] < _cellIndex[0] && boxCoord[1] < _cellIndex[1] &&
                                                boxCoord[2] < _cellIndex[2]) {
                                                adj.push_back(_cellIndex(boxCoord[0], boxCoord[1], boxCoord[2]));
                                            }
                                        }
                                    }

                                    std::sort(adj.begin(), adj.end());
                                    adj.erase(std::unique(std::begin(adj), std::end(adj)), std::end(adj));
                                    adj.erase(std::remove(std::begin(adj), std::end(adj), _cellIndex(i, j, k)),
                                              std::end(adj));

                                    //log::critical("cell {} with n neighbors: {}" ,cellIdx, adj.size());
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
    };

    void update(const util::PerformanceNode &node) {
        auto t = node.timeit();
        setUpBins(node.subnode("setUpBins"));
    }

    virtual void clear() {
        _head.resize(0);
        _list.resize(0);
        _is_set_up = false;
    };

    const util::Index3D &cellIndex() const {
        return _cellIndex;
    };

    const util::Index2D &neighborIndex() const {
        return _cellNeighbors;
    };

    const std::vector<std::size_t> &neighbors() const {
        return _cellNeighborsContent;
    };

    std::size_t *neighborsBegin(std::size_t cellIndex) {
        auto beginIdx = _cellNeighbors(cellIndex, 1_z);
        return &_cellNeighborsContent.at(beginIdx);
    };

    const std::size_t *neighborsBegin(std::size_t cellIndex) const {
        auto beginIdx = _cellNeighbors(cellIndex, 1_z);
        return &_cellNeighborsContent.at(beginIdx);
    };

    std::size_t *neighborsEnd(std::size_t cellIndex) {
        return neighborsBegin(cellIndex) + nNeighbors(cellIndex);
    };

    const std::size_t *neighborsEnd(std::size_t cellIndex) const {
        return neighborsBegin(cellIndex) + nNeighbors(cellIndex);
    };

    std::size_t nNeighbors(std::size_t cellIndex) const {
        return _cellNeighborsContent.at(_cellNeighbors(cellIndex, 0_z));
    };

    data_type &data() {
        return _data.get();
    };

    const data_type &data() const {
        return _data.get();
    };

    scalar maxCutoff() const {
        return _max_cutoff;
    };

    const HEAD &head() const {
        return _head;
    };

    const LIST &list() const {
        return _list;
    };

    template<typename Function>
    void forEachNeighbor(BoxIterator particle, const Function &function) const;

    template<typename Function>
    void forEachNeighbor(BoxIterator particle, std::size_t cell, const Function &function) const;

    bool cellEmpty(std::size_t index) const {
        return _head.at(index) == 0;
    };

    std::size_t cellOfParticle(std::size_t index) const {
        const auto &entry = data().entry_at(index);
        if (entry.deactivated) {
            throw std::invalid_argument("requested deactivated entry");
        }
        const auto &boxSize = _context.get().boxSize();
        if(!(-.5*boxSize[0] <= entry.pos.x && .5*boxSize[0] > entry.pos.x
             && -.5*boxSize[1] <= entry.pos.y && .5*boxSize[1] > entry.pos.y
             && -.5*boxSize[2] <= entry.pos.z && .5*boxSize[2] > entry.pos.z)) {
            throw std::logic_error("CellLinkedList: requested neighbors of particle that was out of bounds.");
        }
        const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / _cellSize.x));
        const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / _cellSize.y));
        const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / _cellSize.z));
        return _cellIndex(i, j, k);
    };

    std::size_t nCells() const {
        return _cellIndex.size();
    };

    BoxIterator particlesBegin(std::size_t cellIndex);

    BoxIterator particlesBegin(std::size_t cellIndex) const;

    BoxIterator particlesEnd(std::size_t cellIndex);

    BoxIterator particlesEnd(std::size_t cellIndex) const;

protected:
    virtual void setUpBins(const util::PerformanceNode &node) {
        if (_max_cutoff > 0) {
            auto t = node.timeit();
            auto tt = node.subnode("allocate").timeit();
            auto nParticles = _data.get().size();
            _head.clear();
            _head.resize(_cellIndex.size());
            _list.resize(0);
            _list.resize(nParticles + 1);
            fillBins(node.subnode("fillBins"));
        }
    };

    void fillBins(const util::PerformanceNode &node) {
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
                _list[pidx] = _head.at(cellIndex);
                _head[cellIndex] = pidx;
            }
            ++pidx;
        }
    }

    HEAD _head;
    // particles, 1-indexed
    LIST _list;


    bool _is_set_up{false};

    scalar _skin{0};
    scalar _max_cutoff{0};
    std::uint8_t _radius;

    Vec3 _cellSize{0, 0, 0};

    util::Index3D _cellIndex;
    // index of size (n_cells x (1 + nAdjacentCells)), where the first element tells how many adj cells are stored
    util::Index2D _cellNeighbors;
    // backing vector of _cellNeighbors index of size (n_cells x (1 + nAdjacentCells))
    std::vector<std::size_t> _cellNeighborsContent;

    std::reference_wrapper<data_type> _data;
    std::reference_wrapper<const readdy::model::Context> _context;
};

class BoxIterator {

    using alloc = std::allocator<std::size_t>;

public:

    using difference_type = typename alloc::difference_type;
    using value_type = typename alloc::value_type;
    using reference = typename alloc::const_reference;
    using pointer = typename alloc::const_pointer;
    using iterator_category = std::forward_iterator_tag;
    using size_type = CellLinkedList::LIST::size_type;

    BoxIterator(const CellLinkedList &ccll, std::size_t state) : _ccll(ccll), _state(state), _val(state - 1) {};

    BoxIterator(const BoxIterator &) = default;

    BoxIterator &operator=(const BoxIterator &) = delete;

    BoxIterator(BoxIterator &&) = default;

    BoxIterator &operator=(BoxIterator &&) = delete;

    ~BoxIterator() = default;

    const BoxIterator operator++(int) {
        BoxIterator tmp(*this);
        operator++();
        return tmp;
    };

    pointer operator->() const {
        return &_val;
    }

    BoxIterator &operator++() {
        _state = _ccll.list().at(_state);
        _val = _state - 1;
        return *this;
    };

    value_type operator*() const {
        return _val;
    };

    bool operator==(const BoxIterator &rhs) const {
        return _state == rhs._state;
    };

    bool operator!=(const BoxIterator &rhs) const {
        return _state != rhs._state;
    };

private:
    const CellLinkedList &_ccll;
    std::size_t _state, _val;
};

inline BoxIterator CellLinkedList::particlesBegin(std::size_t cellIndex) {
    return {*this, _head.at(cellIndex)};
}

inline BoxIterator CellLinkedList::particlesBegin(std::size_t cellIndex) const {
    return {*this, _head.at(cellIndex)};
}

inline BoxIterator CellLinkedList::particlesEnd(std::size_t /*cellIndex*/) const {
    return {*this, 0};
}

inline BoxIterator CellLinkedList::particlesEnd(std::size_t /*cellIndex*/) {
    return {*this, 0};
}

template<typename Function>
inline void CellLinkedList::forEachNeighbor(BoxIterator particle, const Function &function) const {
    forEachNeighbor(particle, cellOfParticle(*particle), function);
}

template<typename Function>
inline void CellLinkedList::forEachNeighbor(BoxIterator particle, std::size_t cell,
                                            const Function &function) const {
    std::for_each(std::next(particle, 1), particlesEnd(cell), [&function, particle](auto x) {
        function(x);
    });
    for (auto itNeighCell = neighborsBegin(cell); itNeighCell != neighborsEnd(cell); ++itNeighCell) {
        std::for_each(particlesBegin(*itNeighCell), particlesEnd(*itNeighCell), function);
    }
}

}
}
}
}
