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
 * @file IOUtils.h
 * @brief << brief description >>
 * @author clonker
 * @date 10.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <h5rd/h5rd.h>

#include <readdy/common/macros.h>

#include "KernelContext.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(ioutils)

struct ReactionInfo {
    const char* name {""};
    std::size_t index {0}; // identify reaction in map of vectors, e.g. for reaction records
    short id {-1}; // global unique reaction id
    std::size_t n_educts {0};
    std::size_t n_products {0};
    scalar rate {0};
    scalar educt_distance {0};
    scalar product_distance {0};
    std::array<particle_type_type, 2> educt_types {{0, 0}};
    std::array<particle_type_type, 2> product_types {{0, 0}};
};
struct ParticleTypeInfo {
    const char* name;
    std::size_t type_id;
    scalar diffusion_constant;
};


std::tuple<h5rd::NativeCompoundType, h5rd::STDCompoundType> getReactionInfoMemoryType(h5rd::Object::ParentFileRef ref);
std::tuple<h5rd::NativeCompoundType, h5rd::STDCompoundType> getParticleTypeInfoType(h5rd::Object::ParentFileRef ref);

void writeSimulationSetup(h5rd::Group &group, const KernelContext &context);
void writeReactionInformation(h5rd::Group &group, const KernelContext &context);
void writeParticleTypeInformation(h5rd::Group &group, const KernelContext &context);

NAMESPACE_END(ioutils)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
