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
#include <ostream>
#include <readdy/common/common.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(_internal)

struct READDY_API ReaDDyVec3 {
public:

    using data_arr = std::array<scalar, 3>;

    union {
        struct { scalar x, y, z; };
        data_arr data;
    };

    ReaDDyVec3() : ReaDDyVec3(0, 0, 0) { };

    ReaDDyVec3(scalar x, scalar y, scalar z) : x(x), y(y), z(z) {};

    explicit ReaDDyVec3(const data_arr &xyz) : data(xyz) {};

    ReaDDyVec3 &operator+=(const ReaDDyVec3 &rhs) {
        for(auto i=0; i < 3; ++i) { data[i] += rhs[i]; }
        return *this;
    };

    ReaDDyVec3 &operator-=(const ReaDDyVec3 &rhs) {
        for(auto i=0; i < 3; ++i) {
            data[i] -= rhs.data[i];
        }
        return *this;
    };

    ReaDDyVec3 &operator*=(scalar a) {
        for(auto i=0; i < 3; ++i) { data[i] *= a; }
        return *this;
    };

    ReaDDyVec3 &operator/=(scalar a) {
        for(auto i=0; i < 3; ++i) { data[i] /= a; }
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
        scalar result {0};
        for(auto i=0; i < 3; ++i) { result += data[i]*data[i]; }
        return result;
    };

    scalar operator[](unsigned int i) const {
        return data[i];
    };

    scalar &operator[](unsigned int i) {
        return data[i];
    };

    ReaDDyVec3& invertElementWise() {
        for(auto i=0; i < 3; ++i) {
            data[i] = static_cast<scalar>(1.) / data[i];
        }
        return *this;
    };

    bool operator==(const ReaDDyVec3 &rhs) const {
        bool result {true};
        for(auto i=0; i < 3; ++i) { result &= data[i] == rhs[i]; }
        return result;
    };

    bool operator!=(const ReaDDyVec3 &rhs) const {
        bool result {true};
        for(auto i=0; i < 3; ++i) { result &= data[i] != rhs[i]; }
        return result;
    };

    friend std::ostream &operator<<(std::ostream &os, const ReaDDyVec3 &vec) {
        os << "ReaDDyVec3(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")";
        return os;
    };

    friend ReaDDyVec3 operator+(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        auto result = ReaDDyVec3(lhs);
        result += rhs;
        return result;
    };

    friend ReaDDyVec3 operator+(const ReaDDyVec3 &lhs, scalar rhs) {
        auto copy = ReaDDyVec3(lhs);
        for(auto i=0U; i < 3; ++i) copy[i] += rhs;
        return copy;
    };

    friend ReaDDyVec3 operator-(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        auto copy = ReaDDyVec3(lhs);
        copy -= rhs;
        return copy;
    };

    friend ReaDDyVec3 operator-(const ReaDDyVec3 &lhs, scalar rhs) {
        return lhs + (-1 * rhs);
    };

    friend ReaDDyVec3 operator/(const ReaDDyVec3 &lhs, scalar rhs) {
        auto copy = ReaDDyVec3(lhs);
        for(auto i=0U; i < 3; ++i) copy[i] /= rhs;
        return copy;
    };

    friend bool operator>=(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        bool result {true};
        for(auto i=0; i < 3; ++i) { result &= lhs[i] >= rhs[i]; }
        return result;
    };

    friend bool operator<=(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        bool result {true};
        for(auto i=0; i < 3; ++i) { result &= lhs[i] <= rhs[i]; }
        return result;
    };

    friend bool operator>(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        return !(lhs <= rhs);
    };

    friend bool operator<(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
        return !(lhs >= rhs);
    };

};

inline readdy::scalar operator*(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs) {
    readdy::scalar result {0};
    for(auto i=0; i < 3; ++i) {
        result += lhs[i]*rhs[i];
    }
    return result;
}

inline ReaDDyVec3 operator*(const ReaDDyVec3 &lhs, const readdy::scalar rhs) {
    return {rhs * lhs[0], rhs * lhs[1], rhs * lhs[2]};
}

inline ReaDDyVec3 operator*(const readdy::scalar rhs, const ReaDDyVec3 &lhs) {
    return lhs * rhs;
}


NAMESPACE_END(_internal)
NAMESPACE_END(readdy)
