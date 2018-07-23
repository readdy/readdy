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
