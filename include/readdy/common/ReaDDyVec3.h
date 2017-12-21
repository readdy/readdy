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

    ReaDDyVec3();

    ReaDDyVec3(scalar x, scalar y, scalar z);

    explicit ReaDDyVec3(const data_arr &xyz);

    ReaDDyVec3 &operator+=(const ReaDDyVec3 &rhs);

    ReaDDyVec3 &operator-=(const ReaDDyVec3 &rhs);

    ReaDDyVec3 &operator*=(scalar a);

    ReaDDyVec3 &operator/=(scalar a);

    ReaDDyVec3 cross(const ReaDDyVec3&) const;

    scalar norm() const;

    scalar normSquared() const;

    scalar operator[](unsigned int i) const;

    scalar &operator[](unsigned int i);

    ReaDDyVec3& invertElementWise();

    bool operator==(const ReaDDyVec3 &rhs) const;

    bool operator!=(const ReaDDyVec3 &rhs) const;

    friend std::ostream &operator<<(std::ostream &, const ReaDDyVec3 &);

    friend ReaDDyVec3 operator+(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs);

    friend ReaDDyVec3 operator+(const ReaDDyVec3 &lhs, scalar rhs);

    friend ReaDDyVec3 operator-(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs);

    friend ReaDDyVec3 operator-(const ReaDDyVec3 &lhs, scalar rhs);

    friend ReaDDyVec3 operator/(const ReaDDyVec3 &lhs, scalar rhs);

    friend bool operator>=(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs);

    friend bool operator<=(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs);

    friend bool operator>(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs);

    friend bool operator<(const ReaDDyVec3 &lhs, const ReaDDyVec3 &rhs);

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
