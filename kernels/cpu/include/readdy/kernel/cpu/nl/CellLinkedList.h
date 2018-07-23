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
 * @file CellLinkedList.h
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright GPL-3
 */

#pragma once

#include <cstddef>
#include <readdy/common/Index.h>
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
                   thread_pool &pool);

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

protected:
    virtual void setUpBins(const util::PerformanceNode &node) = 0;

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
    std::reference_wrapper<thread_pool> _pool;
};

class BoxIterator;

class CompactCellLinkedList : public CellLinkedList {
public:

    using HEAD = std::vector<util::thread::copyable_atomic<std::size_t>>;
    using LIST = std::vector<std::size_t>;
    using entry_cref = const data_type::entry_type &;
    using pair_callback = std::function<void(entry_cref, entry_cref)>;

    using iterator_bounds = std::tuple<std::size_t, std::size_t>;

    CompactCellLinkedList(data_type &data, const readdy::model::Context &context,
                          thread_pool &pool);

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
    void forEachNeighbor(std::size_t particle, std::size_t cell, const Function &function) const;

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

    bool _serial{false};

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

    BoxIterator(const CompactCellLinkedList::LIST &list, std::size_t state)
            : _list(list), _state(state), _val(state - 1) {};

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
        _state = _list.at(_state);
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
    const CompactCellLinkedList::LIST &_list;
    std::size_t _state, _val;
};

inline BoxIterator CompactCellLinkedList::particlesBegin(std::size_t cellIndex) {
    return {_list, (*_head.at(cellIndex)).load()};
}

inline BoxIterator CompactCellLinkedList::particlesBegin(std::size_t cellIndex) const {
    return {_list, (*_head.at(cellIndex)).load()};
}

inline BoxIterator CompactCellLinkedList::particlesEnd(std::size_t /*cellIndex*/) const {
    return {_list, 0};
}

inline BoxIterator CompactCellLinkedList::particlesEnd(std::size_t /*cellIndex*/) {
    return {_list, 0};
}


template<typename Function>
inline void CompactCellLinkedList::forEachNeighbor(std::size_t particle, std::size_t cell,
                                                   const Function &function) const {
    std::for_each(particlesBegin(cell), particlesEnd(cell), [&function, particle](auto x) {
        if (x != particle) function(x);
    });
    for (auto itNeighCell = neighborsBegin(cell); itNeighCell != neighborsEnd(cell); ++itNeighCell) {
        std::for_each(particlesBegin(*itNeighCell), particlesEnd(*itNeighCell), function);
    }
}

}
}
}
}
