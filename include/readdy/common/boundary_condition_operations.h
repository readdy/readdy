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

template<typename Container, typename PBC, int DIM = 3>
inline void fixPosition(Vec3 &vec, const Container &box, const PBC &periodic) {
    for(int d = 0; d < DIM; ++d) {
        if(periodic[d]) {
            while(vec[d] >= .5*box[d]) vec[d] -= box[d];
            while(vec[d] < -.5*box[d]) vec[d] += box[d];
        }
    }
};

template<typename PBC, int DIM=3>
inline Vec3 applyPBC(const Vec3 &in, const scalar *const box, const PBC &periodic) {
    Vec3 out(in);
    fixPosition<const scalar* const, PBC, DIM>(out, box, periodic);
    return out;
};

template<typename Container, typename PBC, int DIM = 3>
inline Vec3 shortestDifference(const Vec3 &lhs, const Vec3 &rhs, const Container &box, const PBC &periodic) {
    auto dv = rhs - lhs;
    for(int i = 0; i < DIM; ++i) {
        if(periodic[i]) {
            if (dv[i] > box[i] * .5) dv[i] -= box[i];
            else if (dv[i] <= -box[i] * .5) dv[i] += box[i];
        }
    }
    return dv;
};

template<typename Container, typename PBC, int DIM=3>
inline scalar distSquared(const Vec3 &lhs, const Vec3 &rhs, const Container &box, const PBC &periodic) {
    auto dv = shortestDifference<Container, PBC, DIM>(lhs, rhs, box, periodic);
    return dv * dv;
}

template<typename Container, typename PBC, int DIM=3>
inline scalar dist(const Vec3 &lhs, const Vec3 &rhs, const Container &box, const PBC &periodic) {
    return std::sqrt(distSquared<Container, PBC, DIM>(lhs, rhs, box, periodic));
}

NAMESPACE_END(bcs)
NAMESPACE_END(readdy)
