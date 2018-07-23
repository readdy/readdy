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
 * A three-dimensional vector.
 *
 * @file ReaDDyVec3.h
 * @brief File containing the definitions for ReaDDyVec3.
 * @author clonker
 * @date 21.04.16
 */

#pragma once

#include <array>
#include <cmath>
#include <string>
#include <algorithm>
#include <type_traits>
#include <sstream>
#include <ostream>
#include "Index.h"
#include "FloatingPoints.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(_internal)

namespace detail {
template<typename T>
using is_arithmetic_type = std::enable_if_t<std::is_arithmetic<T>::value, int>;
}

template<unsigned int N, unsigned int M, typename scalar = double, typename index_type = std::size_t>
class ReaDDyMatrix {
public:
    using data_arr = std::array<scalar, N * M>;
    using size_type = index_type;

    ReaDDyMatrix() {
        std::fill(_data.begin(), _data.end(), 0.);
    }

    explicit ReaDDyMatrix(data_arr data) : _data(data) {}

    ReaDDyMatrix(ReaDDyMatrix &&) noexcept = default;

    ReaDDyMatrix(const ReaDDyMatrix &) = default;

    ReaDDyMatrix &operator=(ReaDDyMatrix &&) noexcept = default;

    ReaDDyMatrix &operator=(const ReaDDyMatrix &) = default;

    ~ReaDDyMatrix() = default;

    const data_arr &data() const {
        return _data;
    }

    data_arr &data() {
        return _data;
    }

    static constexpr std::size_t n() {
        return N;
    }

    static constexpr std::size_t m() {
        return M;
    }

    scalar at(index_type i, index_type j) const {
        return _data.at(_index(i, j));
    }

    scalar &at(index_type i, index_type j) {
        return _data.at(_index(i, j));
    }

    ReaDDyMatrix &operator+=(const ReaDDyMatrix &rhs) {
        std::transform(_data.begin(), _data.end(), rhs._data.begin(), _data.begin(), std::plus<scalar>());
        return *this;
    }

    friend ReaDDyMatrix operator+(ReaDDyMatrix lhs, const ReaDDyMatrix &rhs) {
        lhs += rhs;
        return lhs;
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    ReaDDyMatrix &operator*=(const arithmetic a) {
        std::transform(_data.begin(), _data.end(), _data.begin(), std::bind1st(std::multiplies<arithmetic>(), a));
        return *this;
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend ReaDDyMatrix operator*(ReaDDyMatrix lhs, arithmetic a) {
        lhs *= a;
        return lhs;
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend ReaDDyMatrix operator*(arithmetic a, ReaDDyMatrix rhs) {
        rhs *= a;
        return rhs;
    }

    bool operator==(const ReaDDyMatrix &rhs) const {
        return _data == rhs._data;
    }

    bool operator!=(const ReaDDyMatrix &rhs) const {
        return _data != rhs._data;
    }

    friend std::ostream &operator<<(std::ostream &os, const ReaDDyMatrix &matrix) {
        std::stringstream ss;
        ss << "Matrix " << ReaDDyMatrix::n() << " x " << ReaDDyMatrix::m() << ": " << std::endl;
        for(index_type i = 0; i < ReaDDyMatrix::n(); ++i) {
            for(index_type j = 0; j < ReaDDyMatrix::m(); ++j) {
                ss << matrix.at(i, j) << " ";
            }
            ss << std::endl;
        }
        os << ss.str();
        return os;
    }

private:
    data_arr _data;
    util::Index2D _index{static_cast<std::size_t>(N), static_cast<std::size_t>(M)};
};

template<typename scalar>
using ReaDDyMatrix33 = ReaDDyMatrix<3, 3, scalar, std::size_t>;

template<typename scalar=double>
class READDY_API ReaDDyVec3 {
    static_assert(std::is_arithmetic<scalar>::value, "scalar needs to be an arithmetic type");
public:
    using data_arr = std::array<scalar, 3>;

    union {
        struct {
            scalar x, y, z;
        };
        data_arr data;
    };

    ReaDDyVec3() : ReaDDyVec3(0, 0, 0) {};

    ReaDDyVec3(scalar x, scalar y, scalar z) : x(x), y(y), z(z) {};

    explicit ReaDDyVec3(const data_arr &xyz) : data(xyz) {};

    ReaDDyVec3 &operator+=(const ReaDDyVec3 &rhs) {
        std::transform(data.begin(), data.end(), rhs.data.begin(), data.begin(), std::plus<scalar>());
        return *this;
    };

    ReaDDyVec3 &operator-=(const ReaDDyVec3 &rhs) {
        std::transform(data.begin(), data.end(), rhs.data.begin(), data.begin(), std::minus<scalar>());
        return *this;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    ReaDDyVec3 &operator+=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind1st(std::plus<scalar>(), a));
        return *this;
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    ReaDDyVec3 &operator-=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind2nd(std::minus<scalar>(), a));
        return *this;
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    ReaDDyVec3 &operator*=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind1st(std::multiplies<scalar>(), a));
        return *this;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    ReaDDyVec3 &operator/=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind2nd(std::divides<scalar>(), a));
        return *this;
    };

    ReaDDyVec3 cross(const ReaDDyVec3 &other) const {
        return {
                data[1] * other.data[2] - data[2] * other.data[1],
                data[2] * other.data[0] - data[0] * other.data[2],
                data[0] * other.data[1] - data[1] * other.data[0]
        };
    };

    scalar norm() const {
        return std::sqrt(normSquared());
    };

    scalar normSquared() const {
        return std::inner_product(data.begin(), data.end(), data.begin(), static_cast<scalar>(0));
    };

    scalar operator[](std::size_t i) const {
        return data.at(i);
    };

    scalar &operator[](std::size_t i) {
        return data.at(i);
    };

    ReaDDyVec3 &invertElementWise() {
        for (auto i = 0; i < 3; ++i) {
            data[i] = static_cast<scalar>(1.) / data[i];
        }
        return *this;
    };

    bool operator==(const ReaDDyVec3 &rhs) const {
        return data == rhs.data;
    };

    bool operator!=(const ReaDDyVec3 &rhs) const {
        return data != rhs.data;
    };

    bool almostEquals(const ReaDDyVec3 &rhs) const {
        bool result {true};
        for(std::uint8_t i = 0; i < 3; ++i) {
            const auto fp1 = fp::FloatingPoint<float>(data[i]);
            const auto fp2 = fp::FloatingPoint<float>(rhs.data[i]);
            auto x = fp::FloatingPoint<float>::DistanceBetweenSignAndMagnitudeNumbers(fp1.bits(), fp2.bits());
            result &= fp1.AlmostEquals(fp2);
        }
        return result;
    }

    friend std::ostream &operator<<(std::ostream &os, const ReaDDyVec3 &vec) {
        os << "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")";
        return os;
    };

    friend ReaDDyVec3 operator+(ReaDDyVec3 lhs, const ReaDDyVec3 &rhs) {
        lhs += rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend ReaDDyVec3 operator+(ReaDDyVec3 lhs, arithmetic rhs) {
        lhs += rhs;
        return lhs;
    };

    friend ReaDDyVec3 operator-(ReaDDyVec3 lhs, const ReaDDyVec3 &rhs) {
        lhs -= rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend ReaDDyVec3 operator-(ReaDDyVec3 lhs, arithmetic rhs) {
        lhs -= rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend ReaDDyVec3 operator/(ReaDDyVec3 lhs, arithmetic rhs) {
        lhs /= rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend ReaDDyVec3 operator*(ReaDDyVec3 lhs, const arithmetic rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend bool operator>=(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        bool result{true};
        for (auto i = 0; i < 3; ++i) { result &= lhs[i] >= rhs[i]; }
        return result;
    };

    friend bool operator<=(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        bool result{true};
        for (auto i = 0; i < 3; ++i) { result &= lhs[i] <= rhs[i]; }
        return result;
    };

    friend bool operator>(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        return !(lhs <= rhs);
    };

    friend bool operator<(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        return !(lhs >= rhs);
    };

    friend scalar operator*(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        return std::inner_product(lhs.data.begin(), lhs.data.end(), rhs.data.begin(), static_cast<scalar>(0));
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend ReaDDyVec3 operator*(const arithmetic rhs, ReaDDyVec3 lhs) {
        lhs *= rhs;
        return lhs;
    }

};

NAMESPACE_END(_internal)
NAMESPACE_END(readdy)
