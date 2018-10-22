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
 * @file PotentialRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright BSD-3
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
    using PotentialsO1Collection = std::vector<PotentialOrder1 *>;
    using PotentialsO2Collection = std::vector<PotentialOrder2 *>;

    using PotentialsO1Map = std::unordered_map<ParticleTypeId, PotentialsO1Collection>;
    using PotentialsO2Map = util::particle_type_pair_unordered_map<PotentialsO2Collection>;
    using AltPotentialsO2Map = std::unordered_map<ParticleTypeId, std::unordered_map<ParticleTypeId, PotentialsO2Collection>>;

    explicit PotentialRegistry(std::reference_wrapper<const ParticleTypeRegistry> typeRegistry)
            : _types(typeRegistry) {};

    /**
     * Adds a user defined external potential to the registry.
     * @param potential the potential
     */
    void addUserDefined(potentials::PotentialOrder1 *potential) {
        return _registerO1(potential);
    }

    /**
     * Adds a user defined pair potential to the registry.
     * @param potential the potential
     */
    void addUserDefined(potentials::PotentialOrder2 *potential) {
        return _registerO2(potential);
    }

    /**
     * Register a box potential, which is used to confine particles to a cuboid volume. The energy function
     * increases quadratically with respect to the distance from the cuboid edges, resulting in a
     * harmonic repulsion.
     *
     * @param particleType the particle type for which the box potential should take effect
     * @param forceConstant the force constant determines the strength of repulsion
     * @param origin the coordinate of the lower left corner of the box
     * @param extent the extent from the origin
     */
    void addBox(const std::string &particleType, scalar forceConstant, const Vec3 &origin, const Vec3 &extent) {
        addBox(_types(particleType), forceConstant, origin, extent);
    }
    void addBox(ParticleTypeId particleType, scalar forceConstant, const Vec3 &origin, const Vec3 &extent) {
        auto &pots = _ownPotentialsO1[particleType];
        pots.emplace_back(std::make_shared<Box>(particleType, forceConstant, origin, extent));
        _registerO1(pots.back().get());
    }

    /**
     * Register a harmonic repulsion potential.
     *
     * @param type1 first particle type
     * @param type2 second particle type
     * @param forceConstant the force constant
     * @param interactionDistance the interaction distance
     */
    void addHarmonicRepulsion(const std::string &type1, const std::string &type2, scalar forceConstant,
                              scalar interactionDistance) {
        addHarmonicRepulsion(_types(type1), _types(type2), forceConstant, interactionDistance);
    }
    void addHarmonicRepulsion(ParticleTypeId type1, ParticleTypeId type2, scalar forceConstant,
                              scalar interactionDistance) {
        auto &pots = _ownPotentialsP2[std::tie(type1, type2)];
        pots.emplace_back(std::make_shared<HarmonicRepulsion>(type1, type2, forceConstant, interactionDistance));
        _registerO2(pots.back().get());
    }

    /**
     * Register a weak interaction piecewise harmonic potential.
     * @param particleTypeA particle type A
     * @param particleTypeB particle type B
     * @param forceConstant the force constant
     * @param desiredParticleDistance the distance at which it is most favorable
     *        for the particles to be (w.r.t. this potential)
     * @param depth the depth of the energy well
     * @param noInteractionDistance the distance at which this potential has no effect anymore
     */
    void addWeakInteractionPiecewiseHarmonic(ParticleTypeId type1, ParticleTypeId type2,
                                             scalar forceConstant, scalar desiredDist, scalar depth, scalar cutoff) {
        WeakInteractionPiecewiseHarmonic::Configuration conf{desiredDist, depth, cutoff};
        addWeakInteractionPiecewiseHarmonic(type1, type2, forceConstant, conf);
    }
    void addWeakInteractionPiecewiseHarmonic(const std::string &type1, const std::string &type2,
                                             scalar forceConstant, scalar desiredDist, scalar depth, scalar cutoff) {
        addWeakInteractionPiecewiseHarmonic(_types(type1), _types(type2), forceConstant, desiredDist,
                                            depth, cutoff);
    }
    void
    addWeakInteractionPiecewiseHarmonic(const std::string &type1, const std::string &type2, scalar forceConstant,
                                        const WeakInteractionPiecewiseHarmonic::Configuration &config) {
        addWeakInteractionPiecewiseHarmonic(_types(type1), _types(type2), forceConstant, config);
    }
    void
    addWeakInteractionPiecewiseHarmonic(ParticleTypeId type1, ParticleTypeId type2, scalar forceConstant,
                                        const WeakInteractionPiecewiseHarmonic::Configuration &config) {
        auto &pots = _ownPotentialsP2[std::tie(type1, type2)];
        pots.emplace_back(std::make_shared<WeakInteractionPiecewiseHarmonic>(type1, type2, forceConstant, config));
        _registerO2(pots.back().get());
    }

    /**
    * Constructs a Lennard-Jones-type potential between two particle types A and B (where possibly A = B) of the form
    *
    * \f[ V_{\mbox{LJ}}(r) = k(\epsilon , n, m) \left[ \left(\frac{\sigma}{r}\right)^m - \left(\frac{\sigma}{r}\right)^n \right], \f]
    *
    * where n,m are exponent 1 and 2, respectively, with m > n.
    * If shift == true, it will be defined as
    *
    * \f[ V_{\mbox{LJ, shifted}}(r) = V_{\mbox{LJ}}(r) - V_{\mbox{LJ}}(r_{\mbox{cutoff}}) \f]
    *
    * for r <= cutoffDistance, which makes a difference in energy, but not in force.
    *
    * @param particleType1 particle type A
    * @param particleType2 particle type B
    * @param m first exponent
    * @param n second exponent
    * @param cutoffDistance the cutoff distance
    * @param shift if it should be shifted or not
    * @param epsilon the well depth
    * @param sigma the distance at which the inter-particle potential is zero
    */
    void addLennardJones(const std::string &type1, const std::string &type2, unsigned int m, unsigned int n,
                         scalar cutoff, bool shift, scalar epsilon, scalar sigma) {
        addLennardJones(_types(type1), _types(type2), m, n, cutoff, shift, epsilon, sigma);
    }
    void addLennardJones(ParticleTypeId type1, ParticleTypeId type2, unsigned int m, unsigned int n,
                         scalar cutoff, bool shift, scalar epsilon, scalar sigma) {
        auto &pots = _ownPotentialsP2[std::tie(type1, type2)];
        pots.emplace_back(std::make_shared<LennardJones>(type1, type2, m, n, cutoff, shift, epsilon, sigma));
        _registerO2(pots.back().get());
    }

    /**
     * Constructs a potential that describes screened electrostatics with a hard-core repulsion between two
     * particle types A and B (where possibly A = B) of the form
     *
     * \f[ V(r) = C \frac{\exp(-\kappa r)}{r} + D\left(\frac{\sigma}{r}\right)^n, \f]
     *
     * where the first term is the electrostatic interaction, the constant C has the dimensions of an energy times distance. Its value
     * can be positive or negative and depends on the valencies of the particles (see Debye-Hueckel theory). \f$\kappa\f$ is
     * the inverse screening depth. The second term is a hard-core repulsion, that ensures
     * that the potential does not diverge to negative infinity.
     *
     * @param particleType1 particle type A
     * @param particleType2 particle type B
     * @param electrostaticStrength C
     * @param inverseScreeningDepth \f$\kappa\f$
     * @param repulsionStrength D
     * @param repulsionDistance \f$\sigma\f$
     * @param exponent n
     * @param cutoff the distance from which no energies and forces are calculated further
     */
    void addScreenedElectrostatics(const std::string &particleType1, const std::string &particleType2,
                                   scalar electrostaticStrength, scalar inverseScreeningDepth,
                                   scalar repulsionStrength, scalar repulsionDistance, unsigned int exponent,
                                   scalar cutoff) {
        addScreenedElectrostatics(_types(particleType1), _types(particleType2), electrostaticStrength,
                                  inverseScreeningDepth, repulsionStrength, repulsionDistance, exponent, cutoff);
    }
    void addScreenedElectrostatics(ParticleTypeId particleType1, ParticleTypeId particleType2,
                                   scalar electrostaticStrength, scalar inverseScreeningDepth,
                                   scalar repulsionStrength, scalar repulsionDistance, unsigned int exponent,
                                   scalar cutoff) {
        auto &pots = _ownPotentialsP2[std::tie(particleType1, particleType2)];
        pots.emplace_back(std::make_shared<ScreenedElectrostatics>(particleType1, particleType2, electrostaticStrength,
                                                                   inverseScreeningDepth, repulsionStrength,
                                                                   repulsionDistance, exponent, cutoff));
        _registerO2(pots.back().get());
    }

    /**
     * Register a sphere potential, which is used to confine particles inside or outside a spherical volume.
     * The energy function increases quadratically with respect to the distance from the sphere edge,
     * resulting in a harmonic repulsion.
     *
     * @param particleType the particle type for which the sphere potential should take effect
     * @param forceConstant the force constant determines the strength of repulsion
     * @param origin the center of the sphere
     * @param radius the extent of the sphere
     * @param inclusion if true, the potential will include particles, otherwise exclude them from the volume
     */
    void addSphere(const std::string &particleType, scalar forceConstant, const Vec3 &origin, scalar radius, bool inclusion) {
        addSphere(_types(particleType), forceConstant, origin, radius, inclusion);
    }
    void addSphere(ParticleTypeId particleType, scalar forceConstant, const Vec3 &origin, scalar radius, bool inclusion) {
        auto &pots = _ownPotentialsO1[particleType];
        if (inclusion) {
            pots.emplace_back(std::make_shared<Sphere<true>>(particleType, forceConstant, origin, radius));
        } else {
            pots.emplace_back(std::make_shared<Sphere<false>>(particleType, forceConstant, origin, radius));
        }
        _registerO1(pots.back().get());
    }

    /**
    * Register a spherical barrier potential. For positive height it represents a concentric barrier around the point origin
    * with a certain radius. The potential consists of multiple harmonic snippets.
    *
    * @param particleType the particle type for which the potential should take effect
    * @param origin the center of the sphere
    * @param radius the radius of the sphere
    * @param height the energetic height of the barrier, can be negative
    * @param width width of the barrier, behaves like full-width-at-half-maximum (FWHM)
    */
    void addSphericalBarrier(const std::string &particleType, scalar height, scalar width, const Vec3 &origin,
                             scalar radius) {
        addSphericalBarrier(_types(particleType), height, width, origin, radius);
    }
    void addSphericalBarrier(ParticleTypeId particleType, scalar height, scalar width, const Vec3 &origin,
                             scalar radius) {
        auto &pots = _ownPotentialsO1[particleType];
        pots.emplace_back(std::make_shared<SphericalBarrier>(particleType, height, width, origin, radius));
        _registerO1(pots.back().get());
    }

    /**
     * Register a cylindrical potential that confines particles to its inside or outside using a harmonic potential.
     * The cylinder is defined with an origin (any point on the axis of the cylinder),
     * the normal (unit vector along the axis), and the radius of the cylinder. A boolean flag determines if the
     * potential includes or excludes particles.
     *
     * @param particleType the particle type for which the potential should take effect
     * @param forceConstant the force constant determines the strength of repulsion
     * @param origin any point on the axis of the cylinder
     * @param normal a vector that determines the direction of the cylinder's axis, will be normalized to 1 internally
     * @param radius the radius of the cylinder
     * @param inclusion if true, the potential will include particles, otherwise exclude them from the volume
     */
    void addCylinder(const std::string &particleType, scalar forceConstant, const Vec3 &origin, const Vec3 &normal,
                       scalar radius, bool inclusion) {
        addCylinder(_types(particleType), forceConstant, origin, normal, radius, inclusion);
    }

    void addCylinder(ParticleTypeId particleType, scalar forceConstant, const Vec3 &origin, const Vec3 &normal,
                       scalar radius, bool inclusion) {
        auto &pots = _ownPotentialsO1[particleType];
        if (inclusion) {
            pots.emplace_back(std::make_shared<Cylinder<true>>(particleType, forceConstant, origin, normal, radius));
        } else {
            pots.emplace_back(std::make_shared<Cylinder<false>>(particleType, forceConstant, origin, normal, radius));
        }
        _registerO1(pots.back().get());
    }

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

    void _registerO1(PotentialOrder1 *potential) {
        auto typeId = potential->particleType();
        _potentialsO1[typeId].push_back(potential);
    }

    void _registerO2(PotentialOrder2 *potential) {
        auto type1Id = potential->particleType1();
        auto type2Id = potential->particleType2();
        auto pp = std::tie(type1Id, type2Id);
        _potentialsO2[pp].push_back(potential);
        _alternativeO2Registry[type1Id][type2Id].push_back(potential);
        if (type1Id != type2Id) {
            _alternativeO2Registry[type2Id][type1Id].push_back(potential);
        }
    }

};

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
