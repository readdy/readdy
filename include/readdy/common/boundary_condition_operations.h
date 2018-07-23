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
 * @file boundary_condition_operations.h
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright BSD-3
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
