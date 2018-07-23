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
 * 3xN tensor.
 *
 * @file Vec3Tensor.h
 * @brief Vec3Tensor definition file
 * @author clonker
 * @date 07.02.17
 * @copyright BSD-3
 */

#pragma once

#include <cassert>
#include "readdy/common/macros.h"
#include "ReaDDyVec3.h"
#include "common.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

template<std::size_t N>
class Vec3Tensor {
    static_assert(N > 0, "we need more than one vec3");
public:
    explicit Vec3Tensor(const std::array<Vec3, N>& cols) : cols(cols) {}

    explicit Vec3Tensor(std::array<Vec3, N>&& cols) : cols(std::move(cols)) {}

    Vec3& at(std::size_t i) {
        assert(i < N);
        return cols.at(i);
    }

    const Vec3& at(std::size_t i) const {
        assert(i < N);
        return cols.at(i);
    }

    Vec3& operator[](std::size_t i) {
        assert(i < N);
        return cols[i];
    }

    constexpr Vec3& operator[](std::size_t i) const {
        assert(i < N);
        return cols[i];
    }

    Vec3Tensor<N>& operator+(const Vec3Tensor<N>& rhs) {
        auto it = cols.begin();
        auto rhs_it = rhs.cols.cbegin();
        for(; it != cols.end(); ++it, ++rhs_it) {
            *it += *rhs_it;
        }
        return *this;
    }

    Vec3Tensor<N> operator+(const Vec3Tensor<N>& rhs) const {
        Vec3Tensor<N> copy {*this};
        auto it = copy.cols.begin();
        auto rhs_it = rhs.cols.cbegin();
        for(; it != copy.cols.end(); ++it, ++rhs_it) {
            *it += *rhs_it;
        }
        return copy;
    }

    Vec3Tensor<N>& operator*(const scalar& rhs) {
        for(auto& v : cols) {
            v *= rhs;
        }
        return *this;
    }

    Vec3Tensor<N> operator*(const scalar& rhs) const {
        Vec3Tensor<N> copy {*this};
        copy = copy * rhs;
        return copy;
    }

    template<typename = typename std::enable_if<N == 3>::type >
    Vec3 operator*(const Vec3& rhs) const {
        const auto& v0 = cols.at(0);
        const auto& v1 = cols.at(1);
        const auto& v2 = cols.at(2);
        return {
                v0[0] * rhs[0] + v1[0] * rhs[1] + v2[0] * rhs[2],
                v0[1] * rhs[0] + v1[1] * rhs[1] + v2[1] * rhs[2],
                v0[2] * rhs[0] + v1[2] * rhs[1] + v2[2] * rhs[2]
        };
    }
private:
    std::array<Vec3, N> cols;
};

template<std::size_t N>
Vec3Tensor<N>& operator*(const scalar scalar, Vec3Tensor<N>& tensor) {
    return tensor * scalar;
}

Vec3 operator*(const Vec3& vec, const Vec3Tensor<3>& tensor) {
    const auto& v0 = tensor.at(0);
    const auto& v1 = tensor.at(1);
    const auto& v2 = tensor.at(2);
    return {
            v0 * vec,
            v1 * vec,
            v2 * vec
    };
}

NAMESPACE_END(model)
NAMESPACE_END(readdy)
