/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file boundary_condition_operations.h
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include "common.h"
#include "ReaDDyVec3.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(bcs)



template<bool PX, bool PY, bool PZ>
inline void fixPosition(Vec3 &vec, const scalar dx, const scalar dy, const scalar dz);

template<>
inline void fixPosition<true, true, true>(Vec3 &vec, const scalar dx, const scalar dy, const scalar dz) {
    vec[0] -= floor((vec[0] + .5 * dx) / dx) * dx;
    vec[1] -= floor((vec[1] + .5 * dy) / dy) * dy;
    vec[2] -= floor((vec[2] + .5 * dz) / dz) * dz;
};
template<>
inline void fixPosition<true, true, false>(Vec3 &vec, const scalar dx, const scalar dy, const scalar dz) {
    vec[0] -= floor((vec[0] + .5 * dx) / dx) * dx;
    vec[1] -= floor((vec[1] + .5 * dy) / dy) * dy;
};
template<>
inline void fixPosition<true, false, true>(Vec3 &vec, const scalar dx, const scalar dy, const scalar dz) {
    vec[0] -= floor((vec[0] + .5 * dx) / dx) * dx;
    vec[2] -= floor((vec[2] + .5 * dz) / dz) * dz;
};
template<>
inline void fixPosition<true, false, false>(Vec3 &vec, const scalar dx, const scalar dy, const scalar dz) {
    vec[0] -= floor((vec[0] + .5 * dx) / dx) * dx;
};
template<>
inline void fixPosition<false, true, true>(Vec3 &vec, const scalar dx, const scalar dy, const scalar dz) {
    vec[1] -= floor((vec[1] + .5 * dy) / dy) * dy;
    vec[2] -= floor((vec[2] + .5 * dz) / dz) * dz;
};
template<>
inline void fixPosition<false, true, false>(Vec3 &vec, const scalar dx, const scalar dy, const scalar dz) {
    vec[1] -= floor((vec[1] + .5 * dy) / dy) * dy;
};
template<>
inline void fixPosition<false, false, true>(Vec3 &vec, const scalar dx, const scalar dy, const scalar dz) {
    vec[2] -= floor((vec[2] + .5 * dz) / dz) * dz;
};
template<>
inline void fixPosition<false, false, false>(Vec3 &vec, const scalar dx, const scalar dy, const scalar dz) {};

template<bool PX, bool PY, bool PZ>
inline Vec3 applyPBC(const Vec3 &in, scalar dx, scalar dy, scalar dz);

template<>
inline Vec3
applyPBC<true, false, false>(const Vec3 &in, const scalar dx, const scalar /*unused*/, const scalar /*unused*/) {
    return {in[0] - std::floor((in[0] + static_cast<scalar>(.5) * dx) / dx) * dx, in[1], in[2]};
}

template<>
inline Vec3
applyPBC<false, true, false>(const Vec3 &in, const scalar /*unused*/, const scalar dy, const scalar /*unused*/) {
    return {in[0], in[1] - std::floor((in[1] + static_cast<scalar>(.5) * dy) / dy) * dy, in[2]};
}

template<>
inline Vec3
applyPBC<false, false, true>(const Vec3 &in, const scalar /*unused*/, const scalar /*unused*/, const scalar dz) {
    return {in[0], in[1], in[2] - std::floor((in[2] + static_cast<scalar>(.5) * dz) / dz) * dz};
}

template<>
inline Vec3 applyPBC<true, true, false>(const Vec3 &in, const scalar dx, const scalar dy, const scalar /*unused*/) {
    return {in[0] - std::floor((in[0] + static_cast<scalar>(.5) * dx) / dx) * dx,
            in[1] - std::floor((in[1] + static_cast<scalar>(.5) * dy) / dy) * dy, in[2]};
}

template<>
inline Vec3 applyPBC<true, false, true>(const Vec3 &in, const scalar dx, const scalar /*unused*/, const scalar dz) {
    return {in[0] - std::floor((in[0] + static_cast<scalar>(.5) * dx) / dx) * dx, in[1],
            in[2] - std::floor((in[2] + static_cast<scalar>(.5) * dz) / dz) * dz};
}

template<>
inline Vec3 applyPBC<false, true, true>(const Vec3 &in, const scalar /*unused*/, const scalar dy, const scalar dz) {
    return {in[0], in[1] - std::floor((in[1] + static_cast<scalar>(.5) * dy) / dy) * dy,
            in[2] - std::floor((in[2] + static_cast<scalar>(.5) * dz) / dz) * dz};
}

template<>
inline Vec3 applyPBC<true, true, true>(const Vec3 &in, const scalar dx, const scalar dy, const scalar dz) {
    return {in[0] - std::floor((in[0] + static_cast<scalar>(.5) * dx) / dx) * dx,
            in[1] - std::floor((in[1] + static_cast<scalar>(.5) * dy) / dy) * dy,
            in[2] - std::floor((in[2] + static_cast<scalar>(.5) * dz) / dz) * dz};
}

template<>
inline Vec3 applyPBC<false, false, false>(const Vec3 &in, const scalar /*unused*/,
                                          const scalar /*unused*/, const scalar /*unused*/) {
    return in;
}

template<bool PX, bool PY, bool PZ>
inline Vec3 shortestDifference(const Vec3 &lhs, const Vec3 &rhs, const scalar dx, const scalar dy,
                               const scalar dz);

template<>
inline Vec3 shortestDifference<true, true, true>(const Vec3 &lhs, const Vec3 &rhs, const scalar dx, const scalar dy,
                                                 const scalar dz) {
    auto dv = rhs - lhs;
    if (dv[0] > dx * .5) dv[0] -= dx;
    else if (dv[0] <= -dx * .5) dv[0] += dx;
    if (dv[1] > dy * .5) dv[1] -= dy;
    else if (dv[1] <= -dy * .5) dv[1] += dy;
    if (dv[2] > dz * .5) dv[2] -= dz;
    else if (dv[2] <= -dz * .5) dv[2] += dz;
    return dv;
};
template<>
inline Vec3 shortestDifference<true, true, false>(const Vec3 &lhs, const Vec3 &rhs, const scalar dx, const scalar dy,
                                                  const scalar dz) {
    auto dv = rhs - lhs;
    if (dv[0] > dx * .5) dv[0] -= dx;
    else if (dv[0] <= -dx * .5) dv[0] += dx;
    if (dv[1] > dy * .5) dv[1] -= dy;
    else if (dv[1] <= -dy * .5) dv[1] += dy;
    return dv;
};
template<>
inline Vec3 shortestDifference<true, false, true>(const Vec3 &lhs, const Vec3 &rhs, const scalar dx, const scalar dy,
                                                  const scalar dz) {
    auto dv = rhs - lhs;
    if (dv[0] > dx * .5) dv[0] -= dx;
    else if (dv[0] <= -dx * .5) dv[0] += dx;
    if (dv[2] > dz * .5) dv[2] -= dz;
    else if (dv[2] <= -dz * .5) dv[2] += dz;
    return dv;
};
template<>
inline Vec3 shortestDifference<true, false, false>(const Vec3 &lhs, const Vec3 &rhs, const scalar dx, const scalar dy,
                                                  const scalar dz) {
    auto dv = rhs - lhs;
    if (dv[0] > dx * .5) dv[0] -= dx;
    else if (dv[0] <= -dx * .5) dv[0] += dx;
    return dv;
};
template<>
inline Vec3 shortestDifference<false, true, true>(const Vec3 &lhs, const Vec3 &rhs, const scalar dx, const scalar dy,
                                                  const scalar dz) {
    auto dv = rhs - lhs;
    if (dv[1] > dy * .5) dv[1] -= dy;
    else if (dv[1] <= -dy * .5) dv[1] += dy;
    if (dv[2] > dz * .5) dv[2] -= dz;
    else if (dv[2] <= -dz * .5) dv[2] += dz;
    return dv;
};
template<>
inline Vec3 shortestDifference<false, true, false>(const Vec3 &lhs, const Vec3 &rhs, const scalar dx, const scalar dy,
                                                   const scalar dz) {
    auto dv = rhs - lhs;
    if (dv[1] > dy * .5) dv[1] -= dy;
    else if (dv[1] <= -dy * .5) dv[1] += dy;
    return dv;
};
template<>
inline Vec3 shortestDifference<false, false, true>(const Vec3 &lhs, const Vec3 &rhs, const scalar dx, const scalar dy,
                                                   const scalar dz) {
    auto dv = rhs - lhs;
    if (dv[2] > dz * .5) dv[2] -= dz;
    else if (dv[2] <= -dz * .5) dv[2] += dz;
    return dv;
};
template<>
inline Vec3 shortestDifference<false, false, false>(const Vec3 &lhs, const Vec3 &rhs, const scalar dx, const scalar dy,
                                                    const scalar dz) {
    return rhs - lhs;
};

template<bool PX, bool PY, bool PZ>
inline scalar
distSquared(const Vec3 &lhs, const Vec3 &rhs, const scalar dx, const scalar dy, const scalar dz) {
    auto dv = shortestDifference<PX, PY, PZ>(lhs, rhs, dx, dy, dz);
    return dv * dv;
}


NAMESPACE_END(bcs)
NAMESPACE_END(readdy)