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

#include <readdy/common/macros.h>
#include <readdy/io/Group.h>
#include "KernelContext.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(ioutils)

void writeSimulationSetup(io::Group &group, const KernelContext &context);

struct ReactionInfo {
    const char* name;
    std::size_t index {0}; // identify reaction in map of vectors, e.g. for reaction records
    short id {-1}; // global unique reaction id
    std::size_t n_educts {0};
    std::size_t n_products {0};
    double rate {0};
    double educt_distance {0};
    double product_distance {0};
    std::array<particle_type_type, 2> educt_types {{0, 0}};
    std::array<particle_type_type, 2> product_types {{0, 0}};
};

class ReactionInfoMemoryType : public readdy::io::NativeCompoundType {
    static readdy::io::NativeCompoundType get() {
        static readdy::io::NativeCompoundType type = readdy::io::NativeCompoundTypeBuilder(sizeof(ReactionInfo))
                .insertString("name", offsetof(ReactionInfo, name))
                .insert<decltype(std::declval<ReactionInfo>().index)>("index", offsetof(ReactionInfo, index))
                .insert<decltype(std::declval<ReactionInfo>().id)>("id", offsetof(ReactionInfo, id))
                .insert<decltype(std::declval<ReactionInfo>().n_educts)>("n_educts", offsetof(ReactionInfo, n_educts))
                .insert<decltype(std::declval<ReactionInfo>().n_products)>("n_products", offsetof(ReactionInfo, n_products))
                .insert<decltype(std::declval<ReactionInfo>().rate)>("rate", offsetof(ReactionInfo, rate))
                .insert<decltype(std::declval<ReactionInfo>().educt_distance)>("educt_distance", offsetof(ReactionInfo, educt_distance))
                .insert<decltype(std::declval<ReactionInfo>().product_distance)>("product_distance", offsetof(ReactionInfo, product_distance))
                .insertArray<particle_type_type, 2>("educt_types", offsetof(ReactionInfo, educt_types))
                .insertArray<particle_type_type, 2>("product_types", offsetof(ReactionInfo, product_types))
                .build();
        return type;
    }

public:
    ReactionInfoMemoryType() : NativeCompoundType(get()) {}
};

class ReactionInfoFileType : public readdy::io::STDCompoundType {
public:
    ReactionInfoFileType() : STDCompoundType(ReactionInfoMemoryType()) {}
};

void writeReactionInformation(io::Group &group, const KernelContext &context);

struct ParticleTypeInfo {
    const char* name;
    std::size_t type_id;
    double diffusion_constant;
};

class ParticleTypeInfoMemoryType : public readdy::io::NativeCompoundType {
    static readdy::io::NativeCompoundType get() {
        static readdy::io::NativeCompoundType type = readdy::io::NativeCompoundTypeBuilder(sizeof(ParticleTypeInfo))
                .insertString("name", offsetof(ParticleTypeInfo, name))
                .insert<decltype(std::declval<ParticleTypeInfo>().type_id)>("type_id", offsetof(ParticleTypeInfo, type_id))
                .insert<decltype(std::declval<ParticleTypeInfo>().diffusion_constant)>("diffusion_constant", offsetof(ParticleTypeInfo, diffusion_constant))
                .build();
        return type;
    }

public:
    ParticleTypeInfoMemoryType() : NativeCompoundType(get()) {}
};

class ParticleTypeInfoFileType : public readdy::io::STDCompoundType {
public:
    ParticleTypeInfoFileType() : STDCompoundType(ParticleTypeInfoMemoryType()) {}
};

void writeParticleTypeInformation(io::Group &group, const KernelContext &context);

NAMESPACE_END(ioutils)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
