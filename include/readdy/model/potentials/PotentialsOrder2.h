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

    /*scalar getMaximalForce(scalar kbt) const noexcept override {
        scalar fMax1 = forceConstant * conf.desiredParticleDistance;
        scalar fMax2 = 2 * conf.depthAtDesiredDistance *
                       (conf.noInteractionDistance - conf.desiredParticleDistance);
        return std::max(fMax1, fMax2);
    }*/

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

    scalar getCutoffRadiusSquared() const override {
        return cutoffDistanceSquared;
    }

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

    std::string describe() const override;

    scalar getCutoffRadiusSquared() const override {
        return cutoffSquared;
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
