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
 * @file range.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#ifndef READDY_MAIN_RANGE_H
#define READDY_MAIN_RANGE_H
namespace readdy {
namespace util {
template<typename T>
class range {
public:
    class const_iterator {
        friend class range;

    public:
        T operator*() const { return i; }

        const const_iterator &operator++() {
            ++i;
            return *this;
        }

        const_iterator operator++(int) {
            const_iterator copy(*this);
            ++i;
            return copy;
        }

        bool operator==(const const_iterator &other) const { return i == other.i; }

        bool operator!=(const const_iterator &other) const { return i != other.i; }

    protected:
        const_iterator(T start) : i(start) {}

    private:
        T i;
    };

    const_iterator begin() const { return begin_it; }

    const_iterator end() const { return end_it; }

    range(T begin, T end) : begin_it(begin), end_it(end) {}

private:
    const_iterator begin_it;
    const_iterator end_it;
};
}
}


#endif //READDY_MAIN_RANGE_H
