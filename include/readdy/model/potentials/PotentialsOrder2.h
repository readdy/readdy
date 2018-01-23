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
                      scalar forceConstant, scalar interactionDistance)
            : super(type1, type2), _forceConstant(forceConstant), _interactionDistance(interactionDistance),
              _interactionDistanceSquared(interactionDistance*interactionDistance) {}

    scalar interactionDistance() const {
        return _interactionDistance;
    }

    std::string describe() const override;

    scalar getForceConstant() const {
        return _forceConstant;
    }

    scalar getMaximalForce(scalar kbt) const noexcept override {
        return _forceConstant * getCutoffRadius();
    }

    scalar calculateEnergy(const Vec3 &x_ij) const override {
        auto distanceSquared = x_ij * x_ij;
        if (distanceSquared < _interactionDistanceSquared) {
            distanceSquared = std::sqrt(distanceSquared);
            distanceSquared -= _interactionDistance;
            distanceSquared *= distanceSquared;
            return static_cast<scalar>(0.5) * distanceSquared * getForceConstant();
        }
        return 0;
    }

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override {
        auto squared = x_ij * x_ij;
        if (squared < _interactionDistanceSquared && squared > 0) {
            squared = std::sqrt(squared);
            force += (getForceConstant() * (squared - _interactionDistance)) / squared * x_ij;
        } else {
            // nothing happens
        }
    }

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &x_ij) const override {
        auto squared = x_ij * x_ij;
        if (squared < _interactionDistanceSquared && squared > 0) {
            squared = std::sqrt(squared);
            energy += 0.5 * getForceConstant() * std::pow(squared - _interactionDistance, 2);
            force += (getForceConstant() * (squared - _interactionDistance)) / squared * x_ij;
        } else {
            // nothing happens
        }
    }

    scalar getCutoffRadius() const override {
        return _interactionDistance;
    }

    scalar getCutoffRadiusSquared() const override {
        return _interactionDistanceSquared;
    }

    std::string type() const override;

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

    private:
        friend class WeakInteractionPiecewiseHarmonic;

        const scalar desiredParticleDistance, depthAtDesiredDistance, noInteractionDistance, noInteractionDistanceSquared;
    };

    WeakInteractionPiecewiseHarmonic(particle_type_type type1, particle_type_type type2,
                                     scalar forceConstant, const Configuration &config)
            : super(type1, type2), forceConstant(forceConstant), conf(config) {};

    scalar getMaximalForce(scalar kbt) const noexcept override {
        scalar fMax1 = forceConstant * conf.desiredParticleDistance;
        scalar fMax2 = 2 * conf.depthAtDesiredDistance *
                       (conf.noInteractionDistance - conf.desiredParticleDistance);
        return std::max(fMax1, fMax2);
    }

    scalar calculateEnergy(const Vec3 &x_ij) const override {
        const auto dist = std::sqrt(x_ij * x_ij);
        const auto len_part2 = conf.noInteractionDistance - conf.desiredParticleDistance;
        if (dist < conf.desiredParticleDistance) {
            // repulsive as we are closer than the desired distance
            return static_cast<scalar>(.5) * forceConstant * (dist - conf.desiredParticleDistance) * (dist - conf.desiredParticleDistance) -
                   conf.depthAtDesiredDistance;
        }
        // attractive as we are further (but not too far) apart than the desired distance
        if (dist < conf.desiredParticleDistance + c_::half * len_part2) {
            return c_::half * conf.depthAtDesiredDistance * (c_::one / (c_::half * len_part2)) * (c_::one / (c_::half * len_part2)) *
                   (dist - conf.desiredParticleDistance) * (dist - conf.desiredParticleDistance) -
                   conf.depthAtDesiredDistance;
        }
        // if we are not too far apart but still further than in the previous case, attractive
        if (dist < conf.noInteractionDistance) {
            return -c_::half * conf.depthAtDesiredDistance * (c_::one / (c_::half * len_part2)) * (c_::one / (c_::half * len_part2)) *
                   (dist - conf.noInteractionDistance) * (dist - conf.noInteractionDistance);
        }
        return 0;
    }

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override {
        const auto dist = std::sqrt(x_ij * x_ij);
        const auto len_part2 = conf.noInteractionDistance - conf.desiredParticleDistance;
        scalar  factor = 0;
        if (dist < conf.desiredParticleDistance) {
            // repulsive as we are closer than the desired distance
            factor = -1 * forceConstant * (conf.desiredParticleDistance - dist);
        } else {
            // attractive as we are further (but not too far) apart than the desired distance
            if (dist < conf.desiredParticleDistance + .5 * len_part2) {
                factor = -c_::one * conf.depthAtDesiredDistance * (c_::one / (c_::half * len_part2)) * (c_::one / (c_::half * len_part2)) *
                         (conf.desiredParticleDistance - dist);
            } else {
                // if we are not too far apart but still further than in the previous case, attractive
                if (dist < conf.noInteractionDistance) {
                    factor = conf.depthAtDesiredDistance * (c_::one / (c_::half * len_part2)) * (c_::one / (c_::half * len_part2)) *
                             (conf.noInteractionDistance - dist);
                }
            }
        }
        if (dist > 0 && factor != 0) {
            force += factor * x_ij / dist;
        } else {
            // nothing happens
        }
    }

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &x_ij) const override {
        energy += calculateEnergy(x_ij);
        calculateForce(force, x_ij);
    }

    scalar getCutoffRadius() const override {
        return conf.noInteractionDistance;
    }

    scalar getCutoffRadiusSquared() const override {
        return conf.noInteractionDistanceSquared;
    }

    std::string type() const override;

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

    ~LennardJones() override = default;

    scalar calculateEnergy(const Vec3 &x_ij) const override {
        const auto r = x_ij.norm();
        if (r > cutoffDistance) return 0;
        return shift ? energy(r) - energy(cutoffDistance) : energy(r);
    }

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override {
        const auto norm = x_ij.norm();
        if (norm <= cutoffDistance) {
            force += -1. * k * (1 / (sigma * sigma)) *
                     (m * std::pow(sigma / norm, m + 2) - n * std::pow(sigma / norm, n + 2)) *
                     x_ij;
        }
    }

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &x_ij) const override {
        energy += calculateEnergy(x_ij);
        calculateForce(force, x_ij);
    }

    scalar getCutoffRadius() const override {
        return cutoffDistance;
    }

    scalar getCutoffRadiusSquared() const override {
        return cutoffDistanceSquared;
    }

    scalar getMaximalForce(scalar kbt) const noexcept override {
        return 0;
    };

    std::string type() const override;

protected:
    scalar energy(scalar r) const {
        return k * (std::pow(sigma / r, m) - std::pow(sigma / r, n));
    }

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

    ~ScreenedElectrostatics() override = default;

    scalar calculateEnergy(const Vec3 &x_ij) const override {
        const scalar  distance = x_ij.norm();
        scalar  result = electrostaticStrength * std::exp(-inverseScreeningDepth * distance) / distance;
        result += repulsionStrength * std::pow(repulsionDistance / distance, exponent);
        return result;
    }

    void calculateForce(Vec3 &force, const Vec3 &x_ij) const override {
        auto distance = x_ij.norm();
        auto forceFactor = electrostaticStrength * std::exp(-inverseScreeningDepth * distance);
        forceFactor *= (inverseScreeningDepth / distance + c_::one / std::pow(distance, c_::two));
        forceFactor += repulsionStrength * exponent / repulsionDistance * std::pow( repulsionDistance / distance, exponent + c_::one);
        force += forceFactor * (- c_::one * x_ij / distance);
    }

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &x_ij) const override {
        calculateForce(force, x_ij);
        energy += calculateEnergy(x_ij);
    }

    scalar getCutoffRadius() const override {
        return cutoff;
    }

    std::string describe() const override;

    scalar getCutoffRadiusSquared() const override {
        return cutoffSquared;
    }

    scalar getMaximalForce(scalar kbt) const noexcept override {
        return 0;
    }

    std::string type() const override;

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
