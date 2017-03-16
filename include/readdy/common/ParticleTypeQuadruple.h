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
 * @file ParticleTypeQuadruple.h
 * @brief << brief description >>
 * @author clonker
 * @date 17.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <utility>

#include "macros.h"
#include "common.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)

struct ParticleTypeQuadruple {
    particle_type_type t1, t2, t3, t4;
    ParticleTypeQuadruple(particle_type_type t1, particle_type_type t2, particle_type_type t3, particle_type_type t4) {
        // divide into two sets {a, b} and {c, d} and sort them separately
        if(t1 > t2) {
            std::swap(t1, t2);
        }
        if(t3 > t4) {
            std::swap(t3, t4);
        }
        // the smallest first element of the two subsets is the overall smallest element
        if(t1 > t3) {
            std::swap(t1, t3);
        }
        ParticleTypeQuadruple::t1 = t1;

        // now we have {t1, a}, {b, c} of which we know that {b < a} and {t1 < c}
        // check if b < c or c < b:
        if(t3 > t4) {
            ParticleTypeQuadruple::t2 = t4;
            ParticleTypeQuadruple::t3 = t3;
            ParticleTypeQuadruple::t4 = t2;
        } else {
            ParticleTypeQuadruple::t2 = t3;
            if(t2 > t4) {
                ParticleTypeQuadruple::t3 = t4;
                ParticleTypeQuadruple::t4 = t2;
            } else {
                ParticleTypeQuadruple::t3 = t2;
                ParticleTypeQuadruple::t4 = t4;
            }
        }
    }
};

NAMESPACE_END(util)
NAMESPACE_END(readdy)
