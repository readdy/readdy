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
 * This file contains the index_persistent_vector, a vector structure together with a stack of 'blanks'. Removal of
 * elements will result in a push back onto the stack of their respective indices, rendering them 'blank'. This handling
 * potentially increases the memory requirements but avoids the shift of access indices.
 *
 * @file index_persistent_vector.h
 * @brief Definitions for the index_persistent_vector
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
    /**
     * the backing vector type
     */
    using backing_vector = typename std::vector<T>;
    /**
     * stack of blanks (indices) type
     */
    using blanks = std::vector<std::size_t>;

    /**
     * the size type of this, inherited from the backing vector
     */
    using size_type = typename backing_vector::size_type;
    /**
     * the difference type of this, inherited from the backing vector
     */
    using difference_type = typename backing_vector::difference_type;
    /**
     * the allocator type of this, inherited from the backing vector
     */
    using allocator_type = typename backing_vector::allocator_type;
    /**
     * the value type of this, inherited from the backing vector
     */
    using value_type = typename backing_vector::value_type;

    /**
     * the iterator type, same as backing vector's iterator
     */
    using iterator = typename backing_vector::iterator;
    /**
     * the const iterator type, same as backing vector's const iterator
     */
    using const_iterator = typename backing_vector::const_iterator;
    /**
     * the reverse iterator type, same as backing vector's reverse iterator
     */
    using reverse_iterator = typename backing_vector::reverse_iterator;
    /**
     * the const reverse iterator type, same as backing vector's const reverse iterator
     */
    using const_reverse_iterator = typename backing_vector::const_reverse_iterator;

    /**
     * gives access to the backing vector
     * @return a reference to the backing vector
     */
    backing_vector &data() {
        return _backing_vector;
    };

    /**
     * gives const access to the backing vector
     * @return a const reference to the backing vector
     */
    const backing_vector &data() const {
        return _backing_vector;
    }

    /**
     * the size of this container, including blanks
     * @return the size
     */
    std::size_t size() const {
        return _backing_vector.size();
    }

    /**
     * returns whether the container is empty
     * @return true if it is empty, false otherwise
     */
    bool empty() const {
        return _backing_vector.size() == _blanks.size();
    }

    /**
     * clears this container
     */
    void clear() {
        _backing_vector.clear();
        _blanks.clear();
    }

    /**
     * Performs a push_back. If the blanks stack is empty, the element is simply pushed back to the backing vector,
     * otherwise it is inserted at the index the stack's first element is pointing to, which then is erased.
     * @param val the value to insert
     * @return an iterator pointing to the inserted element
     */
    iterator push_back(T &&val) {
        if (_blanks.empty()) {
            _backing_vector.push_back(std::forward<T>(val));
            return std::prev(_backing_vector.end());
        } else {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _backing_vector.at(idx) = std::move(val);
            return _backing_vector.begin() + idx;
        }
    }

    /**
     * Performs a push_back. Same as the implementation with an r-value.
     * @param val the value to insert
     * @return an iterator pointing to the inserted element
     */
    iterator push_back(const T &val) {
        if (_blanks.empty()) {
            _backing_vector.push_back(val);
            return std::prev(_backing_vector.end());
        } else {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _backing_vector.at(idx) = val;
            return _backing_vector.begin() + idx;
        }
    }

    /**
     * Performs an emplace_back. Works in the same way as push_back.
     * @tparam Args argument types
     * @param args arguments
     * @return an iterator pointed to the emplaced element
     */
    template<typename... Args>
    iterator emplace_back(Args&&... args) {
        if (_blanks.empty()) {
            _backing_vector.emplace_back(std::forward<Args>(args)...);
            return std::prev(_backing_vector.end());
        } else {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _backing_vector.get_allocator().construct(&*_backing_vector.begin() + idx, std::forward<Args>(args)...);
            return _backing_vector.begin() + idx;
        }
    }

    /**
     * Removes an element by deactivating it and pushing the index to the blanks stack.
     * @param pos an iterator pointing to the element that should be erased
     */
    void erase(iterator pos) {
        deactivate(pos);
        _blanks.push_back(pos - _backing_vector.begin());
    }

    /**
     * Erases a range of elements.
     * @param start begin of the range, inclusive
     * @param end end of the range, exclusive
     */
    void erase(iterator start, const_iterator end) {
        auto offset = start - _backing_vector.begin();
        for (auto it = start; it != end; ++it, ++offset) {
            deactivate(it);
            _blanks.push_back(offset);
        }
    }

    /**
     * Yields the number of deactivated elements, i.e., size() - n_deactivated() is the effective size of this
     * container.
     * @return the number of deactivated elements
     */
    blanks::size_type n_deactivated() const {
        return _blanks.size();
    }

    /**
     * The effective size of this container (i.e., the size of the underlying vector minus the number of deactivated
     * elements).
     * @return the effective size
     */
    size_type effective_size() const {
        return size() - n_deactivated();
    }

    /**
     * Yields a reference to the element at the requested index.
     * @param index the index
     * @return a reference to the element
     */
    T &at(size_type index) {
        return _backing_vector.at(index);
    }

    /**
     * Yields a const reference to the element at the requested index.
     * @param index the index
     * @return a const reference to the element
     */
    const T &at(size_type index) const {
        return _backing_vector.at(index);
    }

    /**
     * Yields an iterator pointing to the begin of this container.
     * @return the iterator
     */
    iterator begin() noexcept {
        return _backing_vector.begin();
    }

    /**
     * Yields an iterator pointing to the end of this container.
     * @return the iterator
     */
    iterator end() noexcept {
        return _backing_vector.end();
    }

    /**
     * Yields a const iterator pointing to the begin of this container.
     * @return the iterator
     */
    const_iterator begin() const noexcept {
        return _backing_vector.begin();
    }

    /**
     * Yields a const iterator pointing to the begin of this container.
     * @return the iterator
     */
    const_iterator cbegin() const noexcept {
        return _backing_vector.cbegin();
    }

    /**
     * Yields a const iterator pointing to the end of this container.
     * @return the iterator
     */
    const_iterator end() const noexcept {
        return _backing_vector.end();
    }

    /**
     * Yields a const iterator pointing to the end of this container.
     * @return the iterator
     */
    const_iterator cend() const noexcept {
        return _backing_vector.cend();
    }

    /**
     * Yields a reverse iterator to the begin of the reversed structure.
     * @return the reverse iterator.
     */
    reverse_iterator rbegin() noexcept {
        return _backing_vector.rbegin();
    }

    /**
     * Yields a reverse iterator to the end of the reversed structure.
     * @return the reverse iterator
     */
    reverse_iterator rend() noexcept {
        return _backing_vector.rend();
    }

    /**
     * Yields a const reverse iterator to the begin of the reversed structure.
     * @return the const reverse iterator
     */
    const_reverse_iterator rbegin() const noexcept {
        return _backing_vector.rbegin();
    }

    /**
     * Yields a const reverse iterator to the end of the reversed structure.
     * @return the const reverse iterator
     */
    const_reverse_iterator rend() const noexcept {
        return _backing_vector.rend();
    }

    /**
     * Yields a const reverse iterator to the begin of the reversed structure.
     * @return the const reverse iterator
     */
    const_reverse_iterator crbegin() const noexcept {
        return _backing_vector.crbegin();
    }

    /**
     * Yields a const reverse iterator to the end of the reversed structure.
     * @return the const reverse iterator
     */
    const_reverse_iterator crend() const noexcept {
        return _backing_vector.crend();
    }

private:

    /**
     * Deactivate an element if the backing structure contains raw elements.
     * @tparam Q element type
     * @param it the iterator to the element
     */
    template<typename Q = T>
    typename std::enable_if<detail::can_be_deactivated<Q>::value>::type deactivate(iterator it) {
        it->deactivate();
    }

    /**
     * Deactivate an element if the backing structure contains a pointer type.
     * @tparam Q element type
     * @param it the iterator to the element
     */
    template<typename Q = T>
    typename std::enable_if<detail::can_be_ptr_deactivated<Q>::value>::type deactivate(iterator it) {
        (*it)->deactivate();
    }

    blanks _blanks;
    backing_vector _backing_vector;
};

NAMESPACE_END(util)
NAMESPACE_END(readdy)