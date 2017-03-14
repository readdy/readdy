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
#include "macros.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(numeric)

inline constexpr double pi() { return 3.141592653589793238462643383279502884e+00; }

template<typename T, typename D, typename std::enable_if<std::is_arithmetic<T>::value && std::is_arithmetic<D>::value, int>::type = 0>
inline typename std::make_unsigned<T>::type positive_modulo(T i, D n) {
    using return_t = typename std::make_unsigned<T>::type;
    return static_cast<return_t>((i % n + n) % n);
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

NAMESPACE_END(numeric)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
