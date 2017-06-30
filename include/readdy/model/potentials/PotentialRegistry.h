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
 * @file PotentialRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/macros.h>
#include <readdy/common/ParticleTypeTuple.h>
#include <unordered_set>
#include <readdy/model/ParticleTypeRegistry.h>
#include "PotentialOrder1.h"
#include "PotentialOrder2.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(potentials)

class PotentialRegistry {

    using pot_ptr_vec1 = std::vector<std::unique_ptr<potentials::PotentialOrder1>>;
    using pot_ptr_vec1_external = std::vector<potentials::PotentialOrder1 *>;
    using pot_ptr_vec2 = std::vector<std::unique_ptr<potentials::PotentialOrder2>>;
    using pot_ptr_vec2_external = std::vector<potentials::PotentialOrder2 *>;
public:
    using particle_type_registry_ref = std::reference_wrapper<const ParticleTypeRegistry>;

    PotentialRegistry(particle_type_registry_ref typeRegistry);

    PotentialRegistry(const PotentialRegistry &) = delete;

    PotentialRegistry &operator=(const PotentialRegistry &) = delete;

    PotentialRegistry(PotentialRegistry &&) = delete;

    PotentialRegistry &operator=(PotentialRegistry &&) = delete;

    using potentials_o1 = std::vector<PotentialOrder1 *>;
    using potentials_o2 = std::vector<PotentialOrder2 *>;

    using potential_o1_registry = std::unordered_map<particle_type_type, potentials_o1>;
    using potential_o2_registry = util::particle_type_pair_unordered_map<potentials_o2>;

    const Potential::id_t add_external(potentials::PotentialOrder1 *potential);

    const Potential::id_t add_external(potentials::PotentialOrder2 *potential);

    const Potential::id_t add(std::unique_ptr<PotentialOrder1> potential);

    const Potential::id_t add(std::unique_ptr<PotentialOrder2> potential);

    void remove(const Potential::id_t handle);

    const potentials_o1 &potentials_of(const particle_type_type type) const;

    const potential_o1_registry &potentials_order1() const;

    const potentials_o2 &potentials_of(const particle_type_type t1, const particle_type_type t2) const;

    const potential_o2_registry &potentials_order2() const;

    const potentials_o1 &potentials_of(const std::string &type) const;

    const potentials_o2 &potentials_of(const std::string &t1, const std::string &t2) const;

    void configure();

    void debug_output() const;

private:
    using potential_o1_registry_internal = std::unordered_map<particle_type_type, pot_ptr_vec1>;
    using potential_o2_registry_internal = util::particle_type_pair_unordered_map<pot_ptr_vec2>;

    const ParticleTypeRegistry &typeRegistry;

    potential_o1_registry potentialO1Registry{};
    potential_o2_registry potentialO2Registry{};

    potential_o1_registry_internal potentialO1RegistryInternal{};
    potential_o1_registry potentialO1RegistryExternal{};
    potential_o2_registry_internal potentialO2RegistryInternal{};
    potential_o2_registry potentialO2RegistryExternal{};

    pot_ptr_vec1_external defaultPotentialsO1{};
    pot_ptr_vec2_external defaultPotentialsO2{};

};

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
