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
 * flat iterator impl inspired from https://stackoverflow.com/a/3623597/2871028
 *
 * @file flat_iterator.h
 * @brief iterate over container of containers as if it were just one container
 * @author clonker
 * @date 15.09.17
 * @copyright GNU Lesser General Public License v3.0
 */
#pragma once

#include <iterator>

#include "common.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)

template<typename OuterIt>
class flat_iterator {
public:
    using outer_iterator = OuterIt;
    using inner_iterator = typename OuterIt::value_type::iterator;
    using iterator_category = std::forward_iterator_tag;
    using value_type = typename inner_iterator::value_type;
    using difference_type = typename inner_iterator::difference_type;
    using pointer = typename inner_iterator::pointer;
    using reference = typename inner_iterator::reference;

    flat_iterator() = default;

    explicit flat_iterator(outer_iterator it) : _outerIt(it), _outerEnd(it) {}

    flat_iterator(outer_iterator it, outer_iterator end) : _outerIt(it), _outerEnd(end) {
        if (_outerIt == _outerEnd) return;

        _innerIt = _outerIt->begin();
        advance_past_empty_inner_containers();
    }

    reference operator*() const { return *_innerIt; }

    pointer operator->() const { return &*_innerIt; }

    flat_iterator &operator++() {
        ++_innerIt;
        if (_innerIt == _outerIt->end()) advance_past_empty_inner_containers();
        return *this;
    }

    flat_iterator operator++(int) {
        flat_iterator it(*this);
        ++*this;
        return it;
    }

    friend bool operator==(const flat_iterator &a,
                           const flat_iterator &b) {
        if (a._outerIt != b._outerIt) return false;
        return !(a._outerIt != a._outerEnd &&
                 b._outerIt != b._outerEnd &&
                 a._innerIt != b._innerIt);

    }

    friend bool operator!=(const flat_iterator &a,
                           const flat_iterator &b) {
        return !(a == b);
    }


private:

    void advance_past_empty_inner_containers() {
        while (_outerIt != _outerEnd && _innerIt == _outerEnd->end()) {
            ++_outerIt;
            if (_outerIt != _outerEnd) _innerIt = _outerIt->begin();
        }
    }

    outer_iterator _outerIt;
    outer_iterator _outerEnd;
    inner_iterator _innerIt;
};

template<typename OuterIt>
class flat_const_iterator {
public:
    using outer_iterator = OuterIt;
    using inner_iterator = typename OuterIt::value_type::const_iterator;
    using iterator_category = std::forward_iterator_tag;
    using value_type = typename inner_iterator::value_type;
    using difference_type = typename inner_iterator::difference_type;
    using pointer = typename inner_iterator::pointer;
    using reference = typename inner_iterator::reference;

    flat_const_iterator() = default;

    explicit flat_const_iterator(outer_iterator it) : _outerIt(it), _outerEnd(it) {}

    flat_const_iterator(outer_iterator it, outer_iterator end) : _outerIt(it), _outerEnd(end) {
        if (_outerIt == _outerEnd) return;

        _innerIt = _outerIt->begin();
        advance_past_empty_inner_containers();
    }

    reference operator*() const { return *_innerIt; }

    pointer operator->() const { return &*_innerIt; }

    flat_const_iterator operator+(std::size_t offset) const {
        auto copy = *this;
        copy += offset;
        return copy;
    };

    bool operator<(const flat_const_iterator &rhs) const {
        return _outerIt < rhs._outerIt || (_outerIt == rhs._outerIt && _innerIt < rhs._innerIt);
    };

    bool operator>(const flat_const_iterator &rhs) const {
        return _outerIt > rhs._outerIt || (_outerIt == rhs._outerIt && _innerIt > rhs._innerIt);
    };

    bool operator<=(const flat_const_iterator &rhs) const {
        return !(*this > rhs);
    };

    bool operator>=(const flat_const_iterator &rhs) const {
        return !(*this < rhs);
    }

    flat_const_iterator &operator+=(std::size_t offset) {
        while(offset > std::distance(_innerIt, _outerIt->end()) && _outerIt != _outerEnd) {
            offset -= std::distance(_innerIt, _outerIt->end());
            _innerIt = _outerIt->begin();
            advance_past_empty_inner_containers();
        }
        _innerIt += offset;
        return *this;
    }

    friend flat_const_iterator operator+(std::size_t offset, const flat_const_iterator &it) {
        return it + offset;
    }

    flat_const_iterator &operator++() {
        ++_innerIt;
        if (_innerIt == _outerIt->end()) advance_past_empty_inner_containers();
        return *this;
    }

    flat_const_iterator operator++(int) {
        flat_const_iterator it(*this);
        ++*this;
        return it;
    }

    friend bool operator==(const flat_const_iterator &a,
                           const flat_const_iterator &b) {
        if (a._outerIt != b._outerIt) return false;
        return !(a._outerIt != a._outerEnd &&
                 b._outerIt != b._outerEnd &&
                 a._innerIt != b._innerIt);

    }

    friend bool operator!=(const flat_const_iterator &a,
                           const flat_const_iterator &b) {
        return !(a == b);
    }


private:

    void advance_past_empty_inner_containers() {
        while (_outerIt != _outerEnd && _innerIt == _outerEnd->end()) {
            ++_outerIt;
            if (_outerIt != _outerEnd) _innerIt = _outerIt->begin();
        }
    }

    outer_iterator _outerIt;
    outer_iterator _outerEnd;
    inner_iterator _innerIt;
};

NAMESPACE_END(util)
NAMESPACE_END(readdy)
