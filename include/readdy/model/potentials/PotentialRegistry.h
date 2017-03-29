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

    using rdy_ptp = readdy::util::particle_type_pair;
    using rdy_ptp_hasher = readdy::util::particle_type_pair_hasher;
    using rdy_ptp_eq = readdy::util::particle_type_pair_equal_to;

    using potential_o1_registry_internal = std::unordered_map<particle_type_type, pot_ptr_vec1>;
    using potential_o2_registry_internal = std::unordered_map<rdy_ptp, pot_ptr_vec2, rdy_ptp_hasher, rdy_ptp_eq>;
public:
    using particle_type_registry_ref = std::reference_wrapper<const ParticleTypeRegistry>;

    PotentialRegistry(particle_type_registry_ref typeRegistry);

    using potentials_o1 = std::vector<PotentialOrder1*>;
    using potentials_o2 = std::vector<PotentialOrder2*>;

    using rdy_pot_1_registry = std::unordered_map<particle_type_type, potentials_o1>;

    using rdy_pot_2_registry = std::unordered_map<rdy_ptp, potentials_o2, rdy_ptp_hasher, rdy_ptp_eq>;

    const short registerExternalPotential(potentials::PotentialOrder1 *potential);

    const short registerExternalPotential(potentials::PotentialOrder2 *potential);

    template<typename R>
    const short registerPotential(std::unique_ptr<R> potential,
                                  typename std::enable_if<std::is_base_of<potentials::PotentialOrder1, R>::value>::type * = 0) {
        const auto id = potential->getId();
        auto typeId = typeRegistry.getParticleTypeID(potential->particleType);
        if (potentialO1RegistryInternal.find(typeId) == potentialO1RegistryInternal.end()) {
            potentialO1RegistryInternal.insert(std::make_pair(typeId, pot_ptr_vec1()));
        }
        potentialO1RegistryInternal[typeId].push_back(std::move(potential));
        return id;
    }

    template<typename R>
    const short registerPotential(std::unique_ptr<R> potential,
                                  typename std::enable_if<std::is_base_of<potentials::PotentialOrder2, R>::value>::type * = 0) {
        const auto id = potential->getId();
        auto type1Id = typeRegistry.getParticleTypeID(potential->particleType1);
        auto type2Id = typeRegistry.getParticleTypeID(potential->particleType2);
        auto pp = std::tie(type1Id, type2Id);
        if (potentialO2RegistryInternal.find(pp) == potentialO2RegistryInternal.end()) {
            potentialO2RegistryInternal.emplace(pp, pot_ptr_vec2());
        }
        potentialO2RegistryInternal[pp].push_back(std::move(potential));
        return id;
    }

    void deregisterPotential(const short handle);

    const potentials_o1& getOrder1Potentials(const particle_type_type type) const;

    const rdy_pot_1_registry& getAllOrder1Potentials() const;

    std::vector<particle_type_type> getAllOrder1RegisteredPotentialTypes() const;

    const potentials_o2 &getOrder2Potentials(const particle_type_type t1, const particle_type_type t2) const;

    const rdy_pot_2_registry& getAllOrder2Potentials() const;

    std::vector<util::particle_type_pair> getAllOrder2RegisteredPotentialTypes() const;

    /**
     * Get an unstructured list of all second order potentials. Useful in neighborlists, where maxcutoff is required.
     */
    const potentials_o2 getVectorAllOrder2Potentials() const;

    const potentials_o1& getOrder1Potentials(const std::string &type) const;

    const potentials_o2 &getOrder2Potentials(const std::string &t1, const std::string &t2) const;

    void configure();

    void debug_output() const;

private:
    const ParticleTypeRegistry& typeRegistry;

    rdy_pot_1_registry potentialO1Registry{};
    rdy_pot_2_registry potentialO2Registry{};

    potential_o1_registry_internal potentialO1RegistryInternal{};
    rdy_pot_1_registry potentialO1RegistryExternal{};
    potential_o2_registry_internal potentialO2RegistryInternal{};
    rdy_pot_2_registry potentialO2RegistryExternal{};

    pot_ptr_vec1_external defaultPotentialsO1{};
    pot_ptr_vec2_external defaultPotentialsO2{};

};

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
