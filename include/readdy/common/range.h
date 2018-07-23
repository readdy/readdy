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
