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
 * This header contains definitions for the POD types ReactionInfo and ParticleTypeInfo, which are used to store
 * defined reactions and particle types in a hdf5 file.
 * It also contains getter methods for their respective hdf5 compound types as well as helper methods which will
 * write the simulation setup / the reaction information / the particle type information to file.
 *
 * @file IOUtils.h
 * @brief Definitions for ReactionInfo, ParticleTypeInfo and hdf5 helper methods.
 * @author clonker
 * @date 10.03.17
 * @copyright GPL-3
 */

#pragma once

#include <vector>

#include <h5rd/h5rd.h>

#include <readdy/common/macros.h>
#include "Context.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(ioutils)

struct ReactionInfo {
    const char* name {""};
    // std::size_t index {0}; // identify reaction in map of vectors, e.g. for reaction records
    reactions::Reaction::ReactionId id {0}; // global unique reaction id
    std::size_t n_educts {0};
    std::size_t n_products {0};
    scalar rate {0};
    scalar educt_distance {0};
    scalar product_distance {0};
    std::array<ParticleTypeId, 2> educt_types {{0, 0}};
    std::array<ParticleTypeId, 2> product_types {{0, 0}};
};
struct ParticleTypeInfo {
    const char* name;
    std::size_t type_id;
    scalar diffusion_constant;
};

std::tuple<h5rd::NativeCompoundType, h5rd::STDCompoundType> getReactionInfoMemoryType(h5rd::Object::ParentFileRef ref);
std::tuple<h5rd::NativeCompoundType, h5rd::STDCompoundType> getParticleTypeInfoType(h5rd::Object::ParentFileRef ref);

void writeGeneralContextInformation(h5rd::Group &group, const Context &context);
void writeSimulationSetup(h5rd::Group &group, const Context &context);
void writeReactionInformation(h5rd::Group &group, const Context &context);
void writeParticleTypeInformation(h5rd::Group &group, const Context &context);

NAMESPACE_END(ioutils)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
