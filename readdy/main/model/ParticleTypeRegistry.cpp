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
 * @file ParticleTypeRegistry.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright BSD-3
 */

#include <readdy/model/ParticleTypeRegistry.h>
#include <readdy/model/Utils.h>

#include <utility>

namespace readdy::model {


ParticleTypeInfo::ParticleTypeInfo(std::string name, const DiffusionConstant diffusionConstant,
                                   const ParticleFlavor flavor, const ParticleTypeId typeId)
        : name(std::move(name)), diffusionConstant(diffusionConstant), flavor(flavor), typeId(typeId) {}


void ParticleTypeRegistry::add(const std::string &name, DiffusionConstant diffusionConst,
                               const ParticleFlavor flavor) {
    util::validateTypeName(name);
    {
        if(std::holds_alternative<scalar>(diffusionConst)) {
            if(std::get<scalar>(diffusionConst) < 0) {
                throw std::invalid_argument("The diffusion constant must not be negative");
            }
        } else {
            const auto &darr = std::get<Vec3>(diffusionConst).data;
            if(std::any_of(darr.begin(), darr.end(), [](auto x) {return x < 0;})) {
                throw std::invalid_argument("The diffusion constant must not be negative");
            }
        }
        if (std::holds_alternative<Vec3>(diffusionConst)) {
            auto dvec = std::get<Vec3>(diffusionConst);
            auto fpx = fp::FloatingPoint<scalar>(dvec[0]);
            auto fpy = fp::FloatingPoint<scalar>(dvec[1]);
            auto fpz = fp::FloatingPoint<scalar>(dvec[2]);
            if(fpx.AlmostEquals(fpy) && fpx.AlmostEquals(fpz)) {
                diffusionConst = dvec[0];  // change to single scalar if almost equal
            }
        }
        // check if name already exists
        for(const auto &e : particle_info_) {
            if(e.second.name == name) {
                throw std::invalid_argument(fmt::format("A particle type with name {} already exists.", name));
            }
        }
    }
    ParticleTypeId t_id = type_counter_++;
    type_mapping_.emplace(name, t_id);
    particle_info_.emplace(std::make_pair(t_id, ParticleTypeInfo{name, diffusionConst, flavor, t_id}));
    n_types_++;
}

std::string ParticleTypeRegistry::describe() const {
    namespace rus = readdy::util::str;
    std::string description;
    description += fmt::format(" - particle types:\n");
    auto flavorName = [](auto flavor) -> std::string {
        switch (flavor) {
            case particleflavor::NORMAL: {
                return "";
            }
            case particleflavor::TOPOLOGY: {
                return "Topology";
            }
            case particleflavor::MEMBRANE: {
                return "Membrane";
            }
            default: {
                throw std::logic_error("encountered a particle flavor that was neither of NORMAL, TOPOLOGY, MEMBRANE.");
            }
        }
    };
    for (const auto &entry : particle_info_) {
        if (std::holds_alternative<scalar>(entry.second.diffusionConstant)) {
            description += fmt::format("     * {} particle type \"{}\" with D={}\n", flavorName(entry.second.flavor),
                                       entry.second.name, std::get<scalar>(entry.second.diffusionConstant));
        } else {
            description += fmt::format("     * {} particle type \"{}\" with D={}\n", flavorName(entry.second.flavor),
                                       entry.second.name, std::get<Vec3>(entry.second.diffusionConstant));
        }
    }
    return description;
}

}
