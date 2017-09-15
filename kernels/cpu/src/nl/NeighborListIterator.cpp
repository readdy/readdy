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
 * @file NeighborListIterator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 15.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/cpu/nl/NeighborList.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

const IteratorState::size_type *IteratorState::begin() const {
    if(_adaptive) {
        return &*_adaptiveIt->begin();
    } else {
        return &*(_staticIt+1);
    }
}

const IteratorState::size_type *IteratorState::end() const {
    if(_adaptive) {
        return &*_adaptiveIt->end();
    } else {
        return &*(_staticIt + *_staticIt + 1);
    }
}

bool IteratorState::adaptive() const {
    return _adaptive;
}

bool NeighborListIterator::operator!=(const NeighborListIterator &rhs) const {
    return !(rhs==*this);
}

bool NeighborListIterator::operator==(const NeighborListIterator &rhs) const {
    return _state._adaptive == rhs._state._adaptive &&
            (_state._adaptive ? _state._adaptiveIt == rhs._state._adaptiveIt : _state._staticIt == rhs._state._staticIt);
}

NeighborListIterator::reference NeighborListIterator::operator*() const {
    return _state;
}

NeighborListIterator::pointer NeighborListIterator::operator->() const {
    return &_state;
}

NeighborListIterator &NeighborListIterator::operator++() {
    if(_state._adaptive) {
        ++_state._adaptiveIt;
    } else {
        auto size = *_state._staticIt;
        _state._staticIt += size+1;
    }
    return *this;
}

bool NeighborListIterator::operator<(const NeighborListIterator &rhs) const {
    return _state._adaptive ? _state._adaptiveIt < rhs._state._adaptiveIt : _state._staticIt < rhs._state._staticIt;
}

bool NeighborListIterator::operator>(const NeighborListIterator &rhs) const {
    return _state._adaptive ? _state._adaptiveIt > rhs._state._adaptiveIt : _state._staticIt > rhs._state._staticIt;
}

bool NeighborListIterator::operator<=(const NeighborListIterator &rhs) const {
    return _state._adaptive ? _state._adaptiveIt <= rhs._state._adaptiveIt : _state._staticIt <= rhs._state._staticIt;
}

bool NeighborListIterator::operator>=(const NeighborListIterator &rhs) const {
    return _state._adaptive ? _state._adaptiveIt >= rhs._state._adaptiveIt : _state._staticIt >= rhs._state._staticIt;
}

NeighborListIterator &NeighborListIterator::operator+=(size_type x) {
    if(_state._adaptive) {
        _state._adaptiveIt += x;
    } else {
        for(int i = 0; i < x; ++i) operator++();
    }
    return *this;
}

NeighborListIterator NeighborListIterator::operator+(NeighborListIterator::size_type x) const {
    auto copy = NeighborListIterator(*this);
    copy += x;
    return copy;
}

NeighborListIterator operator+(NeighborListIterator::size_type x, const NeighborListIterator &it) {
    return it + x;
}

NeighborListIterator::NeighborListIterator(IteratorState::static_iterator staticIterator) : _state() {
    _state._adaptive = false;
    _state._staticIt = staticIterator;
}

NeighborListIterator::NeighborListIterator(IteratorState::adaptive_iterator adaptiveIterator) : _state() {
    _state._adaptive = true;
    _state._adaptiveIt = adaptiveIterator;
}

NeighborListIterator::NeighborListIterator(const NeighborListIterator &rhs) = default;

NeighborListIterator &NeighborListIterator::operator=(const NeighborListIterator &) = default;

NeighborListIterator::NeighborListIterator(NeighborListIterator &&) = default;

NeighborListIterator &NeighborListIterator::operator=(NeighborListIterator &&) = default;

}
}
}
}
