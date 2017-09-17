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
 * @file NeighborListIterator.h
 * @brief << brief description >>
 * @author clonker
 * @date 15.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <cstddef>
#include <vector>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class NeighborListIterator;

struct IteratorState {
    using size_type = std::size_t;
    using const_iterator = std::vector<std::vector<std::size_t>>::const_iterator;
    using inner_iterator = std::vector<std::size_t>::const_iterator;

    const size_type *begin() const;
    const size_type *end() const;

    IteratorState() = default;

    IteratorState(const IteratorState &) = default;
    IteratorState(IteratorState &&) = default;

    IteratorState &operator=(const IteratorState &) = default;
    IteratorState &operator=(IteratorState &&) = default;

    ~IteratorState() = default;

    bool adaptive() const;

    std::size_t n_neighbors() const;

    std::size_t current_particle() const;

private:
    friend class NeighborListIterator;

    std::size_t _adaptive_pidx;
    bool _adaptive {true};
    const_iterator _iterator {};
    inner_iterator _inner_iterator {};
};

class NeighborListIterator {
    using Alloc = std::allocator<IteratorState>;
public:
    using difference_type = typename Alloc::difference_type;
    using value_type = typename Alloc::value_type;
    using reference = typename Alloc::const_reference;
    using pointer = typename Alloc::const_pointer;
    using iterator_category = std::forward_iterator_tag;
    using size_type = IteratorState::size_type;

    explicit NeighborListIterator(IteratorState::const_iterator iterator, bool adaptive);

    NeighborListIterator(const NeighborListIterator &);

    NeighborListIterator &operator=(const NeighborListIterator &);

    NeighborListIterator(NeighborListIterator &&) = default;

    NeighborListIterator &operator=(NeighborListIterator &&) = default;

    ~NeighborListIterator() = default;

    bool operator==(const NeighborListIterator &rhs) const;

    bool operator!=(const NeighborListIterator &rhs) const;

    bool operator<(const NeighborListIterator &rhs) const;

    bool operator>(const NeighborListIterator &rhs) const;

    bool operator<=(const NeighborListIterator &rhs) const;

    bool operator>=(const NeighborListIterator &rhs) const;

    NeighborListIterator &operator++();

    NeighborListIterator &operator+=(size_type);

    NeighborListIterator operator+(size_type) const;

    friend NeighborListIterator operator+(size_type, const NeighborListIterator &);

    reference operator*() const;

    pointer operator->() const;

private:
    IteratorState _state;
};

}
}
}
}