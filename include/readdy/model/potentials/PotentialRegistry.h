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
public:
    using PotentialId = Potential::PotentialId;

    using PotentialsO1Collection = std::vector<PotentialOrder1 *>;
    using PotentialsO2Collection = std::vector<PotentialOrder2 *>;

    using PotentialsO1Map = std::unordered_map<ParticleTypeId, PotentialsO1Collection>;
    using PotentialsO2Map = util::particle_type_pair_unordered_map<PotentialsO2Collection>;
    using AltPotentialsO2Map = std::unordered_map<ParticleTypeId, std::unordered_map<ParticleTypeId, PotentialsO2Collection>>;

    explicit PotentialRegistry(std::reference_wrapper<const ParticleTypeRegistry> typeRegistry)
    : _types(typeRegistry) {};

    PotentialRegistry(const PotentialRegistry &) = default;

    PotentialRegistry &operator=(const PotentialRegistry &) = default;

    PotentialRegistry(PotentialRegistry &&) = default;

    PotentialRegistry &operator=(PotentialRegistry &&) = default;

    ~PotentialRegistry() = default;


    PotentialId addUserDefined(potentials::PotentialOrder1 *potential) {
        return _registerO1(potential);
    }

    PotentialId addUserDefined(potentials::PotentialOrder2 *potential) {
        return _registerO2(potential);
    }

    PotentialId addBox(const std::string &particleType, scalar forceConstant, const Vec3 &origin, const Vec3 &extent) {
        return addBox(_types(particleType), forceConstant, origin, extent);
    }

    PotentialId addBox(ParticleTypeId particleType, scalar forceConstant, const Vec3 &origin, const Vec3 &extent) {
        auto &pots = _ownPotentialsO1[particleType];
        pots.emplace_back(std::make_shared<Box>(particleType, forceConstant, origin, extent));
        return _registerO1(pots.back().get());
    }

    PotentialId addHarmonicRepulsion(const std::string &type1, const std::string &type2, scalar forceConstant,
                                 scalar interactionDistance) {
        return addHarmonicRepulsion(_types(type1), _types(type2), forceConstant, interactionDistance);
    }

    PotentialId addHarmonicRepulsion(ParticleTypeId type1, ParticleTypeId type2, scalar forceConstant,
                                 scalar interactionDistance) {
        auto &pots = _ownPotentialsP2[std::tie(type1, type2)];
        pots.emplace_back(std::make_shared<HarmonicRepulsion>(type1, type2, forceConstant, interactionDistance));
        return _registerO2(pots.back().get());
    }

    PotentialId addWeakInteractionPiecewiseHarmonic(ParticleTypeId type1, ParticleTypeId type2,
                                                scalar forceConstant, scalar desiredDist, scalar depth, scalar cutoff) {
        WeakInteractionPiecewiseHarmonic::Configuration conf{desiredDist, depth, cutoff};
        return addWeakInteractionPiecewiseHarmonic(type1, type2, forceConstant, conf);
    }

    PotentialId addWeakInteractionPiecewiseHarmonic(const std::string &type1, const std::string &type2,
                                                scalar forceConstant, scalar desiredDist, scalar depth, scalar cutoff) {
        return addWeakInteractionPiecewiseHarmonic(_types(type1), _types(type2), forceConstant, desiredDist,
                                                   depth, cutoff);
    }

    PotentialId
    addWeakInteractionPiecewiseHarmonic(const std::string &type1, const std::string &type2, scalar forceConstant,
                                        const WeakInteractionPiecewiseHarmonic::Configuration &config) {
        return addWeakInteractionPiecewiseHarmonic(_types(type1), _types(type2), forceConstant, config);
    }

    PotentialId
    addWeakInteractionPiecewiseHarmonic(ParticleTypeId type1, ParticleTypeId type2, scalar forceConstant,
                                        const WeakInteractionPiecewiseHarmonic::Configuration &config) {
        auto &pots = _ownPotentialsP2[std::tie(type1, type2)];
        pots.emplace_back(std::make_shared<WeakInteractionPiecewiseHarmonic>(type1, type2, forceConstant, config));
        return _registerO2(pots.back().get());
    }

    PotentialId addLennardJones(const std::string &type1, const std::string &type2, unsigned int m, unsigned int n,
                            scalar cutoff, bool shift, scalar epsilon, scalar sigma) {
        return addLennardJones(_types(type1), _types(type2), m, n, cutoff, shift, epsilon, sigma);
    }

    PotentialId addLennardJones(ParticleTypeId type1, ParticleTypeId type2, unsigned int m, unsigned int n,
                            scalar cutoff, bool shift, scalar epsilon, scalar sigma) {
        auto &pots = _ownPotentialsP2[std::tie(type1, type2)];
        pots.emplace_back(std::make_shared<LennardJones>(type1, type2, m, n, cutoff, shift, epsilon, sigma));
        return _registerO2(pots.back().get());
    }

    PotentialId addScreenedElectrostatics(const std::string &particleType1, const std::string &particleType2,
                                      scalar electrostaticStrength, scalar inverseScreeningDepth,
                                      scalar repulsionStrength, scalar repulsionDistance, unsigned int exponent,
                                      scalar cutoff) {
        return addScreenedElectrostatics(_types(particleType1), _types(particleType2), electrostaticStrength,
                                         inverseScreeningDepth, repulsionStrength, repulsionDistance, exponent, cutoff);
    }

    PotentialId addScreenedElectrostatics(ParticleTypeId particleType1, ParticleTypeId particleType2,
                                      scalar electrostaticStrength, scalar inverseScreeningDepth,
                                      scalar repulsionStrength, scalar repulsionDistance, unsigned int exponent,
                                      scalar cutoff) {
        auto &pots = _ownPotentialsP2[std::tie(particleType1, particleType2)];
        pots.emplace_back(std::make_shared<ScreenedElectrostatics>(particleType1, particleType2, electrostaticStrength,
                                                                   inverseScreeningDepth, repulsionStrength,
                                                                   repulsionDistance, exponent, cutoff));
        return _registerO2(pots.back().get());
    }

    PotentialId addSphereOut(const std::string &particleType, scalar forceConstant, const Vec3 &origin, scalar radius) {
        return addSphereOut(_types(particleType), forceConstant, origin, radius);
    }

    PotentialId addSphereOut(ParticleTypeId particleType, scalar forceConstant, const Vec3 &origin, scalar radius) {
        auto &pots = _ownPotentialsO1[particleType];
        pots.emplace_back(std::make_shared<SphereOut>(particleType, forceConstant, origin, radius));
        return _registerO1(pots.back().get());
    }

    PotentialId addSphereIn(const std::string &particleType, scalar forceConstant, const Vec3 &origin, scalar radius) {
        return addSphereIn(_types(particleType), forceConstant, origin, radius);
    }

    PotentialId addSphereIn(ParticleTypeId particleType, scalar forceConstant, const Vec3 &origin, scalar radius) {
        auto &pots = _ownPotentialsO1[particleType];
        pots.emplace_back(std::make_shared<SphereIn>(particleType, forceConstant, origin, radius));
        return _registerO1(pots.back().get());
    }

    PotentialId addSphericalBarrier(const std::string &particleType, scalar height, scalar width, const Vec3 &origin,
                                scalar radius) {
        return addSphericalBarrier(_types(particleType), height, width, origin, radius);
    }

    PotentialId addSphericalBarrier(ParticleTypeId particleType, scalar height, scalar width, const Vec3 &origin,
                                scalar radius) {
        auto &pots = _ownPotentialsO1[particleType];
        pots.emplace_back(std::make_shared<SphericalBarrier>(particleType, height, width, origin, radius));
        return _registerO1(pots.back().get());
    }

    void remove(PotentialId handle);

    const PotentialsO1Collection &potentialsOf(const ParticleTypeId type) const {
        static const auto defaultValue = PotentialsO1Collection{};
        auto it = _potentialsO1.find(type);
        return it != std::end(_potentialsO1) ? it->second : defaultValue;
    }

    const PotentialsO1Map &potentialsOrder1() const {
        return _potentialsO1;
    }

    const PotentialsO2Collection &potentialsOf(const ParticleTypeId t1, const ParticleTypeId t2) const {
        static const auto defaultValue = PotentialsO2Collection{};
        auto it = _potentialsO2.find(std::tie(t1, t2));
        return it != std::end(_potentialsO2) ? it->second : defaultValue;
    }

    const AltPotentialsO2Map::value_type::second_type &potentialsOrder2(const ParticleTypeId t) const {
        static const auto defaultValue = AltPotentialsO2Map::value_type::second_type{};
        auto it = _alternativeO2Registry.find(t);
        return it != _alternativeO2Registry.end() ? it->second : defaultValue;
    }

    const PotentialsO2Map &potentialsOrder2() const {
        return _potentialsO2;
    }

    const PotentialsO1Collection &potentialsOf(const std::string &type) const {
        return potentialsOf(_types(type));
    }

    const PotentialsO2Collection &potentialsOf(const std::string &t1, const std::string &t2) const {
        return potentialsOf(_types(t1), _types(t2));
    }

    std::string describe() const;

private:
    using OwnPotentialsO1 = std::vector<std::shared_ptr<potentials::PotentialOrder1>>;
    using OwnPotentialsO2 = std::vector<std::shared_ptr<potentials::PotentialOrder2>>;
    using OwnPotentialsO1Map = std::unordered_map<ParticleTypeId, OwnPotentialsO1>;
    using OwnPotentialsO2Map = util::particle_type_pair_unordered_map<OwnPotentialsO2>;

    std::reference_wrapper<const ParticleTypeRegistry> _types;

    AltPotentialsO2Map _alternativeO2Registry{};
    PotentialsO1Map _potentialsO1{};
    PotentialsO2Map _potentialsO2{};

    OwnPotentialsO1Map _ownPotentialsO1{};
    OwnPotentialsO2Map _ownPotentialsP2{};

    PotentialId _registerO1(PotentialOrder1 *potential) {
        auto typeId = potential->particleType();
        _potentialsO1[typeId].push_back(potential);
        return potential->getId();
    }

    PotentialId _registerO2(PotentialOrder2 *potential) {
        auto id = potential->getId();
        auto type1Id = potential->particleType1();
        auto type2Id = potential->particleType2();
        auto pp = std::tie(type1Id, type2Id);
        _potentialsO2[pp].push_back(potential);
        _alternativeO2Registry[type1Id][type2Id].push_back(potential);
        if(type1Id != type2Id) {
            _alternativeO2Registry[type2Id][type1Id].push_back(potential);
        }
        return id;
    }

};

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)

#include "misc/PotentialRegistry_misc.h"
