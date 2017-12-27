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
 * @file CellLinkedList.h
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <cstddef>
#include <readdy/common/Index.h>
#include <readdy/common/thread/Config.h>
#include <readdy/model/Context.h>
#include <readdy/common/Timer.h>
#include <readdy/common/thread/atomic.h>
#include <readdy/kernel/cpu/data/DefaultDataContainer.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class CellLinkedList {
public:
    using data_type = readdy::kernel::cpu::data::DefaultDataContainer;
    using cell_radius_type = std::uint8_t;

    CellLinkedList(data_type &data, const readdy::model::Context &context,
                   const readdy::util::thread::Config &config);

    void setUp(scalar skin, cell_radius_type radius, const util::PerformanceNode &node);

    virtual void update(const util::PerformanceNode &node) = 0;
    
    virtual void clear() = 0;
    
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

    std::size_t cellOfParticle(std::size_t index) const {
        const auto &entry = data().entry_at(index);
        if(entry.deactivated) {
            throw std::invalid_argument("requested deactivated entry");
        }
        const auto &boxSize = _context.get().boxSize();
        const auto i = static_cast<std::size_t>(std::floor((entry.pos.x + .5 * boxSize[0]) / _cellSize.x));
        const auto j = static_cast<std::size_t>(std::floor((entry.pos.y + .5 * boxSize[1]) / _cellSize.y));
        const auto k = static_cast<std::size_t>(std::floor((entry.pos.z + .5 * boxSize[2]) / _cellSize.z));
        return _cellIndex(i, j, k);
    };

    std::size_t nCells() const {
        return _cellIndex.size();
    };

protected:
    virtual void setUpBins(const util::PerformanceNode &node) = 0;

    bool _is_set_up {false};

    scalar _skin {0};
    scalar _max_cutoff {0};
    std::uint8_t _radius;

    Vec3 _cellSize {0, 0, 0};

    util::Index3D _cellIndex;
    // index of size (n_cells x (1 + nAdjacentCells)), where the first element tells how many adj cells are stored
    util::Index2D _cellNeighbors;
    // backing vector of _cellNeighbors index of size (n_cells x (1 + nAdjacentCells))
    std::vector<std::size_t> _cellNeighborsContent;

    std::reference_wrapper<data_type> _data;
    std::reference_wrapper<const readdy::model::Context> _context;
    std::reference_wrapper<const readdy::util::thread::Config> _config;
};

class BoxIterator;
class MacroBoxIterator;
class NeighborsIterator;

class CompactCellLinkedList : public CellLinkedList {
public:

    using HEAD = std::vector<util::thread::copyable_atomic<std::size_t>>;
    using LIST = std::vector<std::size_t>;
    using entry_cref = const data_type::entry_type&;
    using pair_callback = std::function<void(entry_cref, entry_cref)>;

    using iterator_bounds = std::tuple<std::size_t, std::size_t>;

    CompactCellLinkedList(data_type &data, const readdy::model::Context &context,
                          const util::thread::Config &config);

    void update(const util::PerformanceNode &node) override {
        auto t = node.timeit();
        setUpBins(node.subnode("setUpBins"));
    };

    void clear() override {
        _head.resize(0);
        _list.resize(0);
    };

    BoxIterator particlesBegin(std::size_t cellIndex);

    BoxIterator particlesBegin(std::size_t cellIndex) const;

    BoxIterator particlesEnd(std::size_t cellIndex);

    BoxIterator particlesEnd(std::size_t cellIndex) const;

    MacroBoxIterator macroCellParticlesBegin(std::size_t cellIndex, int skip=-1) const;

    MacroBoxIterator macroCellParticlesEnd(std::size_t cellIndex) const;

    NeighborsIterator cellNeighborsBegin(std::size_t cellIndex) const;

    NeighborsIterator cellNeighborsEnd(std::size_t cellIndex) const;

    const HEAD &head() const {
        return _head;
    };

    const LIST &list() const {
        return _list;
    };

    bool &serial() {
        return _serial;
    };

    const bool &serial() const {
        return _serial;
    };

    template<typename Function>
    void forEachNeighbor(std::size_t particle, const Function &function) const {
        forEachNeighbor(particle, cellOfParticle(particle), function);
    }

    template<typename Function>
    void forEachNeighbor(std::size_t particle, std::size_t cell, const Function& function) const {
        std::for_each(particlesBegin(cell), particlesEnd(cell), [&function, particle](auto x) {
            if(x != particle) function(x);
        });
        for(auto itNeighCell = neighborsBegin(cell); itNeighCell != neighborsEnd(cell); ++itNeighCell) {
            std::for_each(particlesBegin(*itNeighCell), particlesEnd(*itNeighCell), function);
        }
    }

    bool cellEmpty(std::size_t index) const {
        return (*_head.at(index)).load() == 0;
    };
protected:
    void setUpBins(const util::PerformanceNode &node) override;

    template<bool serial>
    void fillBins(const util::PerformanceNode &node);

    HEAD _head;
    // particles, 1-indexed
    LIST _list;

    bool _serial {false};

};

class BoxIterator {

    using alloc = std::allocator<std::size_t>;

public:

    using difference_type = typename alloc::difference_type;
    using value_type = typename alloc::value_type;
    using reference = typename alloc::const_reference;
    using pointer = typename alloc::const_pointer;
    using iterator_category = std::forward_iterator_tag;
    using size_type = CompactCellLinkedList::LIST::size_type;

    BoxIterator(const CompactCellLinkedList &ccll, std::size_t state);

    BoxIterator(const BoxIterator&) = default;
    BoxIterator &operator=(const BoxIterator &) = default;
    BoxIterator(BoxIterator &&) = default;
    BoxIterator &operator=(BoxIterator &&) = default;
    ~BoxIterator() = default;

    BoxIterator &operator++();

    value_type operator*() const;

    bool operator==(const BoxIterator &rhs) const;

    bool operator!=(const BoxIterator &rhs) const;

private:
    std::reference_wrapper<const CompactCellLinkedList> _ccll;
    std::size_t _state;
};

inline BoxIterator CompactCellLinkedList::particlesBegin(std::size_t cellIndex) {
    return {*this, (*_head.at(cellIndex)).load()};
}

inline BoxIterator CompactCellLinkedList::particlesBegin(std::size_t cellIndex) const {
    return {*this, (*_head.at(cellIndex)).load()};
}

inline BoxIterator CompactCellLinkedList::particlesEnd(std::size_t /*cellIndex*/) const {
    return {*this, 0};
}

inline BoxIterator CompactCellLinkedList::particlesEnd(std::size_t /*cellIndex*/) {
    return {*this, 0};
}

class MacroBoxIterator {
    using alloc = std::allocator<std::size_t>;
public:
    using difference_type = typename alloc::difference_type;
    using value_type = typename alloc::value_type;
    using reference = typename alloc::const_reference;
    using pointer = typename alloc::const_pointer;
    using iterator_category = std::forward_iterator_tag;
    using size_type = CompactCellLinkedList::LIST::size_type;

    MacroBoxIterator(const CompactCellLinkedList &ccll, std::size_t centerCell, const std::size_t *currentCell,
                     bool end=false, int skip=-1);

    MacroBoxIterator &operator++();

    reference operator*() const;

    pointer operator->() const;

    bool operator==(const MacroBoxIterator &rhs) const;
    bool operator!=(const MacroBoxIterator &rhs) const;

private:

    const std::size_t *getCurrentCell(const std::size_t *initial) const;

    BoxIterator getCurrentBoxIterator() const;

    const std::size_t *_neighborCellsBegin;
    const std::size_t *_neighborCellsEnd;

    std::reference_wrapper<const CompactCellLinkedList> _ccll;
    std::size_t _centerCell;

    const std::size_t *_currentCell;
    BoxIterator _currentBoxIt;
    BoxIterator _currentBoxEnd;

    value_type _state;

    int _skip;
};

class NeighborsIterator {
    using Alloc = std::allocator<std::size_t>;
public:

    using difference_type = typename Alloc::difference_type;
    using value_type = typename Alloc::value_type;
    using reference = typename Alloc::const_reference;
    using pointer = typename Alloc::const_pointer;
    using iterator_category = std::forward_iterator_tag;
    using size_type = std::size_t;

    /**
     * Iterator that will iterate over all particles of the specified cell times all "neighboring" particles (i.e.,
     * one through its own cell (cellAt==nullptr) and then through the neighboring cells).
     * @param ccll the cell linked list this iterator belongs to
     * @param cell the cell
     * @param cellAt the cell in the "neighboring cells" iteration, nullptr if own cell
     * @param state the state of this cell
     * @param stateNeigh the state of the neighboring iterator
     */
    NeighborsIterator(const CompactCellLinkedList &ccll, std::size_t cell, std::size_t state);

    NeighborsIterator &operator++();

    reference operator*() const;

    reference currentParticle() const;

    bool operator==(const NeighborsIterator &rhs) const;

    bool operator!=(const NeighborsIterator &rhs) const;

    MacroBoxIterator neighborsBegin() const;

    MacroBoxIterator neighborsEnd() const;

private:
    std::size_t _cell;

    value_type _currentValue;

    BoxIterator _innerIterator;

    std::reference_wrapper<const CompactCellLinkedList> _ccll;
};

inline MacroBoxIterator CompactCellLinkedList::macroCellParticlesBegin(std::size_t cellIndex, int skip) const {
    return {*this, cellIndex, nullptr, false, skip};
}

inline MacroBoxIterator CompactCellLinkedList::macroCellParticlesEnd(std::size_t cellIndex) const {
    return {*this, cellIndex, nullptr, true};
}

inline NeighborsIterator CompactCellLinkedList::cellNeighborsBegin(std::size_t cellIndex) const {
    return {*this, cellIndex, (*_head.at(cellIndex)).load()};
}

inline NeighborsIterator CompactCellLinkedList::cellNeighborsEnd(std::size_t cellIndex) const {
    return {*this, cellIndex, 0};
}

}
}
}
}
