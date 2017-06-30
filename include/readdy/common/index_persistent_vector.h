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
 * @file index_persistent_vector.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <vector>
#include <stack>
#include <algorithm>
#include "macros.h"
#include "traits.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)

/*NAMESPACE_BEGIN(detail)

template<typename>
struct sfinae_true : std::true_type{};

template<typename T>
static auto test_deactivate(int) -> sfinae_true<decltype(std::declval<T>().deactivate())>;
template<typename T>
static auto test_deactivate(long) -> std::false_type;

NAMESPACE_END(detail)

template<typename T>
struct can_be_deactivated : decltype(detail::test_deactivate<T>(0)) {};*/

NAMESPACE_BEGIN(detail)
template<typename T, typename = void>
struct can_be_deactivated : std::false_type {};
template<typename T>
struct can_be_deactivated<T, void_t<decltype(std::declval<T>().deactivate())>> : std::true_type {};

template<typename T, typename = void>
struct can_be_ptr_deactivated : std::false_type {};
template<typename T>
struct can_be_ptr_deactivated<T, void_t<decltype(std::declval<T>()->deactivate())>> : std::true_type {};

NAMESPACE_END(detail)

template<typename T>
class index_persistent_vector {
    static_assert(detail::can_be_deactivated<T>::value || detail::can_be_ptr_deactivated<T>::value,
                  "index_persistent_vector can only work with (ptr) element types which have a deactivate() method");
public:
    using backing_vector = typename std::vector<T>;
    using blanks = std::vector<std::size_t>;

    using size_type = typename backing_vector::size_type;
    using difference_type = typename backing_vector::difference_type;
    using allocator_type = typename backing_vector::allocator_type;
    using value_type = typename backing_vector::value_type;

    using iterator = typename backing_vector::iterator;
    using const_iterator = typename backing_vector::const_iterator;
    using reverse_iterator = typename backing_vector::reverse_iterator;
    using const_reverse_iterator = typename backing_vector::const_reverse_iterator;

    backing_vector &data() {
        return _backing_vector;
    };

    const backing_vector &data() const {
        return _backing_vector;
    }

    std::size_t size() const {
        return _backing_vector.size();
    }

    bool empty() const {
        return _backing_vector.size() == _blanks.size();
    }

    void clear() {
        _backing_vector.clear();
        _blanks.clear();
    }

    iterator push_back(T &&val) {
        if (_blanks.empty()) {
            _backing_vector.push_back(std::forward<T>(val));
            return _backing_vector.end()-1;
        } else {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _backing_vector.at(idx) = std::move(val);
            return _backing_vector.begin() + idx;
        }
    }

    iterator push_back(const T &val) {
        if (_blanks.empty()) {
            _backing_vector.push_back(val);
            return _backing_vector.end()-1;
        } else {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _backing_vector.at(idx) = val;
            return _backing_vector.begin() + idx;
        }
    }

    void erase(iterator pos) {
        deactivate(pos);
        _blanks.push_back(pos - _backing_vector.begin());
    }

    void erase(iterator start, const_iterator end) {
        auto offset = start - _backing_vector.begin();
        for (auto it = start; it != end; ++it, ++offset) {
            deactivate(it);
            _blanks.push_back(offset);
        }
    }

    blanks::size_type n_deactivated() const {
        return _blanks.size();
    }

    T &at(size_type index) {
        return _backing_vector.at(index);
    }

    const T &at(size_type index) const {
        return _backing_vector.at(index);
    }

    iterator begin() noexcept {
        return _backing_vector.begin();
    }

    iterator end() noexcept {
        return _backing_vector.end();
    }

    const_iterator begin() const noexcept {
        return _backing_vector.begin();
    }

    const_iterator cbegin() const noexcept {
        return _backing_vector.cbegin();
    }

    const_iterator end() const noexcept {
        return _backing_vector.end();
    }

    const_iterator cend() const noexcept {
        return _backing_vector.cend();
    }

    reverse_iterator rbegin() noexcept {
        return _backing_vector.rbegin();
    }

    reverse_iterator rend() noexcept {
        return _backing_vector.rend();
    }

    const_reverse_iterator rbegin() const noexcept {
        return _backing_vector.rbegin();
    }

    const_reverse_iterator rend() const noexcept {
        return _backing_vector.rend();
    }

    const_reverse_iterator crbegin() const noexcept {
        return _backing_vector.crbegin();
    }

    const_reverse_iterator crend() const noexcept {
        return _backing_vector.crend();
    }

private:

    template<typename Q = T>
    typename std::enable_if<detail::can_be_deactivated<Q>::value>::type deactivate(iterator it) {
        it->deactivate();
    }

    template<typename Q = T>
    typename std::enable_if<detail::can_be_ptr_deactivated<Q>::value>::type deactivate(iterator it) {
        (*it)->deactivate();
    }

    blanks _blanks;
    backing_vector _backing_vector;
};

NAMESPACE_END(util)
NAMESPACE_END(readdy)