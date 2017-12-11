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

#include <readdy/kernel/cpu_legacy/nl/NeighborList.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

const IteratorState::size_type *IteratorState::begin() const {
    if(_adaptive) {
        return &*_iterator->begin();
    } else {
        return &*(_inner_iterator+2);
    }
}

const IteratorState::size_type *IteratorState::end() const {
    if(_adaptive) {
        return &*_iterator->end();
    } else {
        return &*(_inner_iterator + n_neighbors() + 2);
    }
}

bool IteratorState::adaptive() const {
    return _adaptive;
}

std::size_t IteratorState::n_neighbors() const {
    return adaptive() ? _iterator->size() : *(_inner_iterator+1);
}

std::size_t IteratorState::current_particle() const {
    if (adaptive()) {
        return _adaptive_pidx;
    } else {
        return *_inner_iterator;
    }
}

bool NeighborListIterator::operator!=(const NeighborListIterator &rhs) const {
    return !(rhs==*this);
}

bool NeighborListIterator::operator==(const NeighborListIterator &rhs) const {
    if(rhs._state._adaptive == _state._adaptive) {
        if(_state._adaptive) {
            return _state._iterator == rhs._state._iterator;
        } else {
            if(_state._iterator == rhs._state._iterator) {
                return _state._iterator == _globalEnd || _state._inner_iterator == rhs._state._inner_iterator;
            }
            return false;
        }
    }
    return false;
}

NeighborListIterator::reference NeighborListIterator::operator*() const {
    return _state;
}

NeighborListIterator::pointer NeighborListIterator::operator->() const {
    return &_state;
}

NeighborListIterator &NeighborListIterator::operator++() {
    if(_state._adaptive) {
        ++_state._iterator;
        ++_state._adaptive_pidx;
    } else {
        _state._inner_iterator = _state._inner_iterator + _state.n_neighbors() + 2;
        while( _state._iterator != _globalEnd && _state._inner_iterator == _state._iterator->end()) {
            ++_state._iterator;
            if(_state._iterator != _globalEnd) {
                _state._inner_iterator = _state._iterator->begin();
            }
        }
    }
    return *this;
}

bool NeighborListIterator::operator<(const NeighborListIterator &rhs) const {
    return _state._iterator < rhs._state._iterator && _state._adaptive ? true : _state._inner_iterator < rhs._state._inner_iterator;
}

bool NeighborListIterator::operator>(const NeighborListIterator &rhs) const {
    return _state._iterator < rhs._state._iterator && _state._adaptive ? true : _state._inner_iterator > rhs._state._inner_iterator;
}

bool NeighborListIterator::operator<=(const NeighborListIterator &rhs) const {
    return !(*this > rhs);
}

bool NeighborListIterator::operator>=(const NeighborListIterator &rhs) const {
    return !(*this < rhs);
}

NeighborListIterator &NeighborListIterator::operator+=(size_type x) {
    if(_state._adaptive) {
        _state._iterator += x;
        _state._adaptive_pidx += x;
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

NeighborListIterator::NeighborListIterator(input_iterator iteratorBegin, input_iterator iteratorEnd, bool adaptive)
        : _state(), _globalEnd(iteratorEnd) {
    _state._adaptive = adaptive;
    _state._iterator = iteratorBegin;
    if(!adaptive) {
        while(_state._iterator != _globalEnd && _state._iterator->begin() == _state._iterator->end()) {
            ++_state._iterator;
        }
        if(_state._iterator != _globalEnd) {
            _state._inner_iterator = _state._iterator->begin();
        }
    } else {
        _state._adaptive_pidx = 0;
    }
}

NeighborListIterator::NeighborListIterator(const NeighborListIterator &rhs) = default;

NeighborListIterator &NeighborListIterator::operator=(const NeighborListIterator &) = default;

}
}
}
}
