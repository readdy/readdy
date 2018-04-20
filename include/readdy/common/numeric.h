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
 * @file numeric.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.09.16
 */

#pragma once

#include <type_traits>

namespace readdy {
namespace util {
namespace numeric {

template<typename scalar>
inline constexpr scalar pi() { return 3.141592653589793238462643383279502884e+00; }

template<typename T, typename D, typename std::enable_if<
        std::is_arithmetic<T>::value && std::is_arithmetic<D>::value, int>::type = 0>
inline typename std::make_unsigned<T>::type positive_modulo(T i, D n) {
    using return_type = typename std::make_unsigned<T>::type;
    return static_cast<return_type>((i % n + n) % n);
}

template<typename T, typename D, typename std::enable_if<
        std::is_arithmetic<T>::value && std::is_arithmetic<D>::value, int>::type = 0>
inline T signed_positive_modulo(T i, D n) {
    return static_cast<T>((i % n + n) % n);
}

template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
constexpr T clamp(T d, T min, T max) {
    const T t = d < min ? min : d;
    return t > max ? max : t;
}

template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
constexpr T clamp_min(T d, T min) {
    return d < min ? min : d;
}
}
}
namespace math {
template<typename Matrix33, typename Vec3>
inline Matrix33 outerProduct(const Vec3 &lhs, const Vec3 &rhs) {
    typename Matrix33::data_arr result;
    for(typename Matrix33::size_type i = 0; i < Matrix33::n(); ++i) {
        for(typename Matrix33::size_type j = 0; j < Matrix33::m(); ++j) {
            result.at(Matrix33::m()*j+i) = lhs[j]*rhs[i];
        }
    }
    return Matrix33(result);
}
}
}
