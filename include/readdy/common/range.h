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
 * header file containing the readdy::util::range() class, which can be used to generate index sets
 *
 * @file range.h
 * @brief range header file
 * @author clonker
 * @date 20.10.16
 */

#pragma once
NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
/**
 * the range class, which can be used to create ranges
 * @tparam T the type of the range
 */
template<typename T>
class range {
public:
    /**
     * a corresponding const_iterator so that iteration through the range is possible
     */
    class const_iterator {
        // we are friends with range
        friend class range;

    public:
        /**
         * get the current element
         * @return the current element
         */
        T operator*() const { return i; }

        /**
         * advance the iterator by one
         * @return myself
         */
        const const_iterator &operator++() {
            ++i;
            return *this;
        }

        /**
         * advance the iterator by one
         * @return myself
         */
        const_iterator operator++(int) {
            const_iterator copy(*this);
            ++i;
            return copy;
        }

        /**
         * equality operator
         * @param other the other iterator
         * @return true if the other iterator points to the same element of the range
         */
        bool operator==(const const_iterator &other) const { return i == other.i; }

        /**
         * neq operator
         * @param other the other iterator
         * @return true if the other iterator does not point to the same element
         */
        bool operator!=(const const_iterator &other) const { return i != other.i; }

    protected:
        /**
         * constructs a new iterator
         * @param start where to start
         */
        const_iterator(T start) : i(start) {}

    private:
        T i;
    };

    /**
     * constructs a new const_iterator at the begin of the range
     * @return the const iterator
     */
    const_iterator begin() const { return begin_it; }

    /**
     * constructs a new const_iterator at the end of the range
     * @return the const iterator
     */
    const_iterator end() const { return end_it; }

    /**
     * constructs a new range [include, end)
     * @param begin the begin (inclusive)
     * @param end the end (exclusive)
     */
    range(T begin, T end) : begin_it(begin), end_it(end) {}

private:
    const_iterator begin_it;
    const_iterator end_it;
};
NAMESPACE_END(util)
NAMESPACE_END(readdy)
