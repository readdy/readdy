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
#include "NeighborListContainer.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class NeighborListIterator;

struct IteratorState {
    using size_type = std::size_t;
    using static_iterator = NeighborListContainer::const_iterator;
    using adaptive_iterator = std::vector<std::vector<std::size_t>>::const_iterator;

    const size_type *begin() const;
    const size_type *end() const;

    IteratorState() = default;

    IteratorState(const IteratorState &) = default;
    IteratorState(IteratorState &&) = default;

    IteratorState &operator=(const IteratorState &) = default;
    IteratorState &operator=(IteratorState &&) = default;

    ~IteratorState() = default;

    bool adaptive() const;

private:
    friend class NeighborListIterator;

    bool _adaptive {true};
    static_iterator _staticIt;
    adaptive_iterator _adaptiveIt;
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

    explicit NeighborListIterator(IteratorState::static_iterator staticIterator);

    explicit NeighborListIterator(IteratorState::adaptive_iterator adaptiveIterator);

    NeighborListIterator(const NeighborListIterator &);

    NeighborListIterator &operator=(const NeighborListIterator &);

    NeighborListIterator(NeighborListIterator &&);

    NeighborListIterator &operator=(NeighborListIterator &&);

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