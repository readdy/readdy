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
 * @file Vec3.h
 * @brief File containing the definitions for Vec3.
 * @author clonker
 * @date 21.04.16
 */

#ifndef READDY_MAIN_VEC3_H
#define READDY_MAIN_VEC3_H

#include <array>
#include <math.h>
#include <readdy/common/macros.h>

namespace readdy {
namespace model {

class READDY_API Vec3 {
public:

    using entry_t = double;

    Vec3();

    Vec3(entry_t x, entry_t y, entry_t z);

    Vec3(const std::array<entry_t, 3> &xyz);

    Vec3 &operator+=(const Vec3 &rhs);

    Vec3 &operator*=(const entry_t a);

    Vec3 &operator/=(const entry_t a);

    entry_t operator[](const unsigned int i) const;

    entry_t &operator[](const unsigned int i);

    bool operator==(const Vec3 &rhs) const;

    bool operator!=(const Vec3 &rhs) const;

    friend std::ostream &operator<<(std::ostream &, const Vec3 &);

    friend Vec3 operator+(const Vec3 &lhs, const Vec3 &rhs);

    friend Vec3 operator+(const Vec3 &lhs, const entry_t rhs);

    friend Vec3 operator-(const Vec3 &lhs, const Vec3 &rhs);

    friend Vec3 operator-(const Vec3 &lhs, const entry_t rhs);

    friend Vec3 operator/(const Vec3 &lhs, const entry_t rhs);

    friend bool operator>=(const Vec3 &lhs, const Vec3 &rhs);

    friend bool operator<=(const Vec3 &lhs, const Vec3 &rhs);

    friend bool operator>(const Vec3 &lhs, const Vec3 &rhs);

    friend bool operator<(const Vec3 &lhs, const Vec3 &rhs);

private:
    std::array<entry_t, 3> data;
};

inline Vec3::entry_t operator*(const Vec3 &lhs, const Vec3 &rhs) {
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

inline Vec3 operator*(const Vec3 &lhs, const Vec3::entry_t rhs) {
    return Vec3(rhs * lhs[0], rhs * lhs[1], rhs * lhs[2]);
}

inline Vec3 operator*(const Vec3::entry_t rhs, const Vec3 &lhs) {
    return lhs * rhs;
}

template<bool PX, bool PY, bool PZ>
inline void fixPosition(Vec3 &vec, const Vec3::entry_t dx, const Vec3::entry_t dy, const Vec3::entry_t dz) {
    if (PX) {
        vec[0] -= floor((vec[0] + .5 * dx) / dx) * dx;
    }
    if (PY) {
        vec[1] -= floor((vec[1] + .5 * dy) / dy) * dy;
    }
    if (PZ) {
        vec[2] -= floor((vec[2] + .5 * dz) / dz) * dz;
    }
};

template<bool PX, bool PY, bool PZ>
inline Vec3 shortestDifference(const Vec3 &lhs, const Vec3 &rhs, const Vec3::entry_t dx, const Vec3::entry_t dy,
                               const Vec3::entry_t dz) {
    auto dv = rhs - lhs;
    if (PX) {
        if (dv[0] > dx * .5) dv[0] -= dx;
        else if (dv[0] <= -dx * .5) dv[0] += dx;
    }
    if (PY) {
        if (dv[1] > dy * .5) dv[1] -= dy;
        else if (dv[1] <= -dy * .5) dv[1] += dy;
    }
    if (PZ) {
        if (dv[2] > dz * .5) dv[2] -= dz;
        else if (dv[2] <= -dz * .5) dv[2] += dz;
    }
    return dv;
};

template<bool PX, bool PY, bool PZ>
inline Vec3::entry_t
distSquared(const Vec3 &lhs, const Vec3 &rhs, const Vec3::entry_t dx, const Vec3::entry_t dy, const Vec3::entry_t dz) {
    auto dv = shortestDifference<PX, PY, PZ>(lhs, rhs, dx, dy, dz);
    return dv * dv;
};

}
}

#endif //READDY_MAIN_VEC3_H
