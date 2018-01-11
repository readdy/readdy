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
#include <unordered_set>

#include <readdy/common/ParticleTypeTuple.h>

#include <readdy/model/ParticleTypeRegistry.h>

#include "PotentialOrder1.h"
#include "PotentialOrder2.h"
#include "PotentialsOrder2.h"
#include "PotentialsOrder1.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(potentials)

class PotentialRegistry {

    using pot_ptr_vec1 = std::vector<std::shared_ptr<potentials::PotentialOrder1>>;
    using pot_ptr_vec1_external = std::vector<potentials::PotentialOrder1 *>;
    using pot_ptr_vec2 = std::vector<std::shared_ptr<potentials::PotentialOrder2>>;
    using pot_ptr_vec2_external = std::vector<potentials::PotentialOrder2 *>;
public:
    using particle_type_registry_ref = std::reference_wrapper<const ParticleTypeRegistry>;

    using id_type = Potential::id_type;

    explicit PotentialRegistry(particle_type_registry_ref typeRegistry) : _types(typeRegistry) {};

    PotentialRegistry(const PotentialRegistry &) = default;

    PotentialRegistry &operator=(const PotentialRegistry &) = default;

    PotentialRegistry(PotentialRegistry &&) = default;

    PotentialRegistry &operator=(PotentialRegistry &&) = default;

    ~PotentialRegistry() = default;

    using potentials_o1 = std::vector<PotentialOrder1 *>;
    using potentials_o2 = std::vector<PotentialOrder2 *>;

    using potential_o1_registry = std::unordered_map<particle_type_type, potentials_o1>;
    using potential_o2_registry = util::particle_type_pair_unordered_map<potentials_o2>;
    using o2_registry_alt = std::unordered_map<particle_type_type, std::unordered_map<particle_type_type, potentials_o2>>;

    id_type addUserDefined(potentials::PotentialOrder1 *potential);

    id_type addUserDefined(potentials::PotentialOrder2 *potential);

    id_type addBox(const std::string &particleType, scalar forceConstant, const Vec3 &origin, const Vec3 &extent) {
        return addBox(_types.get()(particleType), forceConstant, origin, extent);
    }

    id_type addBox(particle_type_type particleType, scalar forceConstant, const Vec3 &origin, const Vec3 &extent) {
        auto &pots = potentialO1RegistryInternal[particleType];
        pots.emplace_back(std::make_shared<Box>(particleType, forceConstant, origin, extent));
        return pots.back()->getId();
    }

    id_type addHarmonicRepulsion(const std::string &type1, const std::string &type2, scalar forceConstant,
                                 scalar interactionDistance) {
        return addHarmonicRepulsion(_types.get()(type1), _types.get()(type2), forceConstant, interactionDistance);
    }

    id_type addHarmonicRepulsion(particle_type_type type1, particle_type_type type2, scalar forceConstant,
                                 scalar interactionDistance) {
        auto &pots = potentialO2RegistryInternal[std::tie(type1, type2)];
        pots.emplace_back(std::make_shared<HarmonicRepulsion>(type1, type2, forceConstant, interactionDistance));
        return pots.back()->getId();
    }

    id_type addWeakInteractionPiecewiseHarmonic(particle_type_type type1, particle_type_type type2,
                                                scalar forceConstant, scalar desiredDist, scalar depth, scalar cutoff) {
        WeakInteractionPiecewiseHarmonic::Configuration conf{desiredDist, depth, cutoff};
        return addWeakInteractionPiecewiseHarmonic(type1, type2, forceConstant, conf);
    }

    id_type addWeakInteractionPiecewiseHarmonic(const std::string &type1, const std::string &type2,
                                                scalar forceConstant, scalar desiredDist, scalar depth, scalar cutoff) {
        return addWeakInteractionPiecewiseHarmonic(_types.get()(type1), _types.get()(type2), forceConstant, desiredDist,
                                                   depth, cutoff);
    }

    id_type
    addWeakInteractionPiecewiseHarmonic(const std::string &type1, const std::string &type2, scalar forceConstant,
                                        const WeakInteractionPiecewiseHarmonic::Configuration &config) {
        return addWeakInteractionPiecewiseHarmonic(_types.get()(type1), _types.get()(type2), forceConstant, config);
    }

    id_type
    addWeakInteractionPiecewiseHarmonic(particle_type_type type1, particle_type_type type2, scalar forceConstant,
                                        const WeakInteractionPiecewiseHarmonic::Configuration &config) {
        auto &pots = potentialO2RegistryInternal[std::tie(type1, type2)];
        pots.emplace_back(std::make_shared<WeakInteractionPiecewiseHarmonic>(type1, type2, forceConstant, config));
        return pots.back()->getId();
    }

    id_type addLennardJones(const std::string &type1, const std::string &type2, unsigned int m, unsigned int n,
                            scalar cutoff, bool shift, scalar epsilon, scalar sigma) {
        return addLennardJones(_types.get()(type1), _types.get()(type2), m, n, cutoff, shift, epsilon, sigma);
    }

    id_type addLennardJones(particle_type_type type1, particle_type_type type2, unsigned int m, unsigned int n,
                            scalar cutoff, bool shift, scalar epsilon, scalar sigma) {
        auto &pots = potentialO2RegistryInternal[std::tie(type1, type2)];
        pots.emplace_back(std::make_shared<LennardJones>(type1, type2, m, n, cutoff, shift, epsilon, sigma));
        return pots.back()->getId();
    }

    id_type addScreenedElectrostatics(const std::string &particleType1, const std::string &particleType2,
                                      scalar electrostaticStrength, scalar inverseScreeningDepth,
                                      scalar repulsionStrength, scalar repulsionDistance, unsigned int exponent,
                                      scalar cutoff) {
        return addScreenedElectrostatics(_types.get()(particleType1), _types.get()(particleType2), electrostaticStrength,
                                         inverseScreeningDepth, repulsionStrength, repulsionDistance, exponent, cutoff);
    }

    id_type addScreenedElectrostatics(particle_type_type particleType1, particle_type_type particleType2,
                                      scalar electrostaticStrength, scalar inverseScreeningDepth,
                                      scalar repulsionStrength, scalar repulsionDistance, unsigned int exponent,
                                      scalar cutoff) {
        auto &pots = potentialO2RegistryInternal[std::tie(particleType1, particleType2)];
        pots.emplace_back(std::make_shared<ScreenedElectrostatics>(particleType1, particleType2, electrostaticStrength,
                                                                   inverseScreeningDepth, repulsionStrength,
                                                                   repulsionDistance, exponent, cutoff));
        return pots.back()->getId();
    }

    id_type addSphereOut(const std::string &particleType, scalar forceConstant, const Vec3 &origin, scalar radius) {
        return addSphereOut(_types.get()(particleType), forceConstant, origin, radius);
    }

    id_type addSphereOut(particle_type_type particleType, scalar forceConstant, const Vec3 &origin, scalar radius) {
        auto &pots = potentialO1RegistryInternal[particleType];
        pots.emplace_back(std::make_shared<SphereOut>(particleType, forceConstant, origin, radius));
        return pots.back()->getId();
    }

    id_type addSphereIn(const std::string &particleType, scalar forceConstant, const Vec3 &origin, scalar radius) {
        return addSphereIn(_types.get()(particleType), forceConstant, origin, radius);
    }

    id_type addSphereIn(particle_type_type particleType, scalar forceConstant, const Vec3 &origin, scalar radius) {
        auto &pots = potentialO1RegistryInternal[particleType];
        pots.emplace_back(std::make_shared<SphereIn>(particleType, forceConstant, origin, radius));
        return pots.back()->getId();
    }

    id_type addSphericalBarrier(const std::string &particleType, scalar height, scalar width, const Vec3 &origin,
                                scalar radius) {
        return addSphericalBarrier(_types.get()(particleType), height, width, origin, radius);
    }

    id_type addSphericalBarrier(particle_type_type particleType, scalar height, scalar width, const Vec3 &origin,
                                scalar radius) {
        auto &pots = potentialO1RegistryInternal[particleType];
        pots.emplace_back(std::make_shared<SphericalBarrier>(particleType, height, width, origin, radius));
        return pots.back()->getId();
    }

    void remove(Potential::id_type handle);

    const potentials_o1 &potentialsOf(const particle_type_type type) const {
        auto it = potentialO1Registry.find(type);
        return it != potentialO1Registry.end() ? it->second : defaultPotentialsO1;
    }

    const potential_o1_registry &potentialsOrder1() const {
        return potentialO1Registry;
    }

    const potentials_o2 &potentialsOf(const particle_type_type t1, const particle_type_type t2) const {
        auto it = potentialO2Registry.find(std::tie(t1, t2));
        return it != potentialO2Registry.end() ? it->second : defaultPotentialsO2;
    }

    const o2_registry_alt::value_type::second_type &potentialsOrder2(const particle_type_type t) const {
        auto it = _alternativeO2Registry.find(t);
        return it != _alternativeO2Registry.end() ? it->second : defaultAlt;
    }

    const potential_o2_registry &potentialsOrder2() const {
        return potentialO2Registry;
    }

    const potentials_o1 &potentialsOf(const std::string &type) const {
        return potentialsOf(_types.get().idOf(type));
    }

    const potentials_o2 &potentialsOf(const std::string &t1, const std::string &t2) const {
        return potentialsOf(_types.get().idOf(t1), _types.get().idOf(t2));
    }

    void configure();

    std::string describe() const;

private:
    using potential_o1_registry_internal = std::unordered_map<particle_type_type, pot_ptr_vec1>;
    using potential_o2_registry_internal = util::particle_type_pair_unordered_map<pot_ptr_vec2>;

    std::reference_wrapper<const ParticleTypeRegistry> _types;

    potential_o1_registry potentialO1Registry{};
    potential_o2_registry potentialO2Registry{};
    o2_registry_alt _alternativeO2Registry{};

    potential_o1_registry_internal potentialO1RegistryInternal{};
    potential_o1_registry potentialO1RegistryExternal{};
    potential_o2_registry_internal potentialO2RegistryInternal{};
    potential_o2_registry potentialO2RegistryExternal{};

    pot_ptr_vec1_external defaultPotentialsO1{};
    pot_ptr_vec2_external defaultPotentialsO2{};

    o2_registry_alt::value_type::second_type defaultAlt{};

};

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)

#include "misc/PotentialRegistry_misc.h"
