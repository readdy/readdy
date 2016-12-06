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
 * make_unique of https://isocpp.org/files/papers/N3656.txt.
 *
 * @file make_unique.h
 * @brief reference implementation of make_unique, see N3656
 * @author clonker
 * @date 29.04.16
 */

#ifndef READDY_MAIN_MAKE_UNIQUE_H_H
#define READDY_MAIN_MAKE_UNIQUE_H_H

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

namespace std {
template<class T>
struct _Unique_if {
    typedef unique_ptr <T> _Single_object;
};

template<class T>
struct _Unique_if<T[]> {
    typedef unique_ptr<T[]> _Unknown_bound;
};

template<class T, size_t N>
struct _Unique_if<T[N]> {
    typedef void _Known_bound;
};

template<class T, class... Args>
typename _Unique_if<T>::_Single_object
make_unique(Args &&... args) {
    return unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<class T>
typename _Unique_if<T>::_Unknown_bound
make_unique(size_t n) {
    typedef typename remove_extent<T>::type U;
    return unique_ptr<T>(new U[n]());
}

template<class T, class... Args>
typename _Unique_if<T>::_Known_bound
make_unique(Args &&...) = delete;
}

#endif //READDY_MAIN_MAKE_UNIQUE_H_H
