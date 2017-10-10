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
 * This header contains the declarations of order 2 potentials. Currently:
 *   - Harmonic repulsion
 *   - Weak interaction piecewise harmonic
 *
 * @file PotentialsOrder2.h
 * @brief Contains the declaration of order 2 potentials.
 * @author clonker
 * @date 09.06.16
 */

#pragma once

#include <ostream>
#include "PotentialOrder2.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(potentials)

class HarmonicRepulsion : public PotentialOrder2 {
    using super = PotentialOrder2;
public:
    HarmonicRepulsion(particle_type_type type1, particle_type_type type2,
                      scalar forceConstant, scalar interactionDistance);

    scalar interactionDistance() const;

    std::string describe() const override;

    scalar getForceConstant() const;

    scalar getMaximalForce(scalar kbt) const noexcept override;

    scalar calculateEnergy(const Vec3 &x_ij) const override;

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override;

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &x_ij) const override;

    scalar getCutoffRadius() const override;

    scalar getCutoffRadiusSquared() const override;

protected:
    scalar _interactionDistance;
    scalar _interactionDistanceSquared;
    scalar _forceConstant;
};

class WeakInteractionPiecewiseHarmonic : public PotentialOrder2 {
    using super = PotentialOrder2;
public:
    std::string describe() const override;

    class Configuration {
    public:
        Configuration(scalar desiredParticleDistance, scalar depthAtDesiredDistance, scalar noInteractionDistance);

        friend std::ostream &operator<<(std::ostream &os, const Configuration &configuration);

    private:
        friend class WeakInteractionPiecewiseHarmonic;

        const scalar desiredParticleDistance, depthAtDesiredDistance, noInteractionDistance, noInteractionDistanceSquared;
    };

    WeakInteractionPiecewiseHarmonic(particle_type_type type1, particle_type_type type2,
                                     scalar forceConstant, const Configuration &config);

    scalar getMaximalForce(scalar kbt) const noexcept override;

    scalar calculateEnergy(const Vec3 &x_ij) const override;

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override;

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &x_ij) const override;

    scalar getCutoffRadius() const override;

    scalar getCutoffRadiusSquared() const override;

protected:
    const Configuration conf;
    const scalar forceConstant;
};

/**
 * Lennard-Jones potential class
 */
class LennardJones : public PotentialOrder2 {
    using super = PotentialOrder2;
public:
    /**
     * Constructs a Lennard-Jones-type potential between two particle types A and B (where possibly A = B) of the
     * form
     *
     * V_LJ(r) = k(epsilon, n, m) [ (sigma/r)^m - (sigma/r)^n ],
     *
     * where n,m are exponent 1 and 2, respectively, with m > n.
     * If shift == true, it will be defined as
     *
     * V_LJ_shifted(r) = V_LJ(r) - V_LJ(cutoffDistance)
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
    LennardJones(particle_type_type type1, particle_type_type type2,
                 unsigned int m, unsigned int n, scalar cutoffDistance,
                 bool shift, scalar epsilon, scalar sigma);

    std::string describe() const override;

    LennardJones(const LennardJones &) = default;

    LennardJones &operator=(const LennardJones &) = delete;

    LennardJones(LennardJones &&) = default;

    LennardJones &operator=(LennardJones &&) = delete;

    ~LennardJones() override;

    scalar calculateEnergy(const Vec3 &x_ij) const override;

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override;

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &x_ij) const override;

    scalar getCutoffRadius() const override;

    scalar getCutoffRadiusSquared() const override;

    scalar getMaximalForce(scalar kbt) const noexcept override;

protected:
    scalar energy(scalar r) const;

    scalar m, n;
    scalar cutoffDistance, cutoffDistanceSquared;
    bool shift; // V_LJ_trunc = V_LJ(r) - V_LJ(r_cutoff)
    scalar epsilon; // depth
    scalar sigma;
    scalar k;
};

class ScreenedElectrostatics : public PotentialOrder2 {
    using super = PotentialOrder2;
public:
    ScreenedElectrostatics(particle_type_type type1, particle_type_type type2, scalar electrostaticStrength,
                           scalar inverseScreeningDepth, scalar repulsionStrength, scalar repulsionDistance,
                           unsigned int exponent, scalar cutoff);

    ScreenedElectrostatics(const ScreenedElectrostatics &) = default;

    ScreenedElectrostatics &operator=(const ScreenedElectrostatics &) = delete;

    ScreenedElectrostatics(ScreenedElectrostatics &&) = delete;

    ScreenedElectrostatics &operator=(ScreenedElectrostatics &&) = delete;

    ~ScreenedElectrostatics() override;

    scalar calculateEnergy(const Vec3 &x_ij) const override;

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override;

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &x_ij) const override;

    scalar getCutoffRadius() const override;

    std::string describe() const override;

    scalar getCutoffRadiusSquared() const override;

    scalar getMaximalForce(scalar kbt) const noexcept override;

protected:
    scalar electrostaticStrength;
    scalar inverseScreeningDepth;
    scalar repulsionStrength;
    scalar repulsionDistance;
    scalar exponent;
    scalar cutoff;
    scalar cutoffSquared;
};

template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<HarmonicRepulsion, T>::value>::type * = 0) {
    return "HarmonicRepulsion";
}

template<typename T>
const std::string
getPotentialName(typename std::enable_if<std::is_base_of<WeakInteractionPiecewiseHarmonic, T>::value>::type * = 0) {
    return "WeakInteractionPiecewiseHarmonic";
}

template<typename T>
const std::string
getPotentialName(typename std::enable_if<std::is_base_of<LennardJones, T>::value>::type * = 0) {
    return "LennardJones";
}

template<typename T>
const std::string
getPotentialName(typename std::enable_if<std::is_base_of<ScreenedElectrostatics, T>::value>::type * = 0) {
    return "ScreenedElectrostatics";
}
NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
