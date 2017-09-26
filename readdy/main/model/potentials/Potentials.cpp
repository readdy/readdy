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
 * @file Potentials.cpp
 * @brief Core library implementation of order 1 potentials
 * @author clonker
 * @author chrisfroe
 * @date 20.06.16
 */

#include <readdy/model/Kernel.h>
#include <readdy/model/potentials/PotentialsOrder1.h>

namespace readdy {
namespace model {
namespace potentials {

short Potential::counter = 0;

/////////////////////////////////////////////////////////////////////////////
//
// Potentials order 1
//
/////////////////////////////////////////////////////////////////////////////

/**
 * Box potential
 */

Vec3 getMinExtent(const Vec3 &origin, const Vec3 &extent) {
    Vec3 result{0, 0, 0};
    for (auto i = 0; i < 3; i++) {
        if (extent[i] > 0) {
            result[i] = origin[i];
        } else {
            result[i] = origin[i] + extent[i];
        }
    }
    return result;
}

Vec3 getMaxExtent(const Vec3 &origin, const Vec3 &extent) {
    Vec3 result{0, 0, 0};
    for (auto i = 0; i < 3; i++) {
        if (extent[i] > 0) {
            result[i] = origin[i] + extent[i];
        } else {
            result[i] = origin[i];
        }
    }
    return result;
}

Box::Box(particle_type_type particleType, scalar forceConstant, const Vec3 &origin,
                             const Vec3 &extent)
        : super(particleType), origin(origin), extent(extent), forceConstant(forceConstant),
          min(getMinExtent(origin, extent)), max(getMaxExtent(origin, extent)) {}

const Vec3 &Box::getOrigin() const { return origin; }

const Vec3 &Box::getExtent() const { return extent; }

scalar Box::getForceConstant() const { return forceConstant; }

scalar Box::getMaximalForce(scalar /*ignored*/) const noexcept {
    return 0;
}

scalar Box::getRelevantLengthScale() const noexcept {
    return std::min(extent[0], std::min(extent[1], extent[2]));
}

scalar Box::calculateEnergy(const Vec3 &position) const {
    scalar energy = 0;

    for (auto i = 0; i < 3; ++i) {
        if (position[i] < min[i] || position[i] > max[i]) {
            if (position[i] < min[i]) {
                energy += 0.5 * forceConstant * (position[i] - min[i]) * (position[i] - min[i]);
            } else {
                energy += 0.5 * forceConstant * (position[i] - max[i]) * (position[i] - max[i]);
            }
        }
    }

    return energy;
}

void Box::calculateForce(Vec3 &force, const Vec3 &position) const {
    for (auto i = 0; i < 3; i++) {
        if (position[i] < min[i] || position[i] > max[i]) {
            if (position[i] < min[i]) {
                force[i] += -1 * forceConstant * (position[i] - min[i]);
            } else {
                force[i] += -1 * forceConstant * (position[i] - max[i]);
            }
        }
    }
}

void Box::calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const {
    energy += calculateEnergy(position);
    calculateForce(force, position);
}

std::string Box::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<Box>() << "[type: " << _particleType << ", origin: " << origin
       << ", extent: " << extent << ", min: " << min << ", max: " << max << ", forceConstant: "
       << forceConstant << "]";
    return ss.str();
}

/*
 * Sphere Potentials
 */

scalar SphereIn::getRelevantLengthScale() const noexcept {
    return radius;
}

scalar SphereIn::getMaximalForce(scalar /*ignored*/) const noexcept {
    return 0;
}

SphereIn::SphereIn(particle_type_type particleType, scalar f, const Vec3 &origin, scalar radius)
        : super(particleType), origin(origin), radius(radius), forceConstant(f) {}

scalar SphereIn::calculateEnergy(const Vec3 &position) const {
    auto difference = position - origin;
    scalar distanceFromOrigin = difference.norm();
    scalar distanceFromSphere = distanceFromOrigin - radius;
    auto energy = static_cast<scalar>(0.);
    if (distanceFromSphere > 0) {
        energy = static_cast<scalar>(0.5 * forceConstant * distanceFromSphere * distanceFromSphere);
    }
    return energy;
}

void SphereIn::calculateForce(Vec3 &force, const Vec3 &position) const {
    auto difference = position - origin;
    scalar distanceFromOrigin = difference.norm();
    scalar distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere > 0) {
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}

void SphereIn::calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const {
    auto difference = position - origin;
    scalar distanceFromOrigin = difference.norm();
    scalar distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere > 0) {
        energy += 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}

std::string SphereIn::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<SphereIn>() << "[type: " << _particleType << ", origin: " << origin
       << ", radius: " << radius << ", forceConstant: " << forceConstant << "]";
    return ss.str();
}

SphereOut::SphereOut(particle_type_type particleType, scalar forceConstant, const Vec3 &origin, scalar radius)
        : super(particleType), forceConstant(forceConstant), origin(origin), radius(radius) {}

std::string SphereOut::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<SphereIn>() << "[type: " << _particleType << ", origin: " << origin
       << ", radius: " << radius << ", forceConstant: " << forceConstant << "]";
    return ss.str();
}

scalar SphereOut::getRelevantLengthScale() const noexcept  {
    return radius;
}

scalar SphereOut::getMaximalForce(scalar /*kbt*/) const noexcept {
    return 0;
}

scalar SphereOut::calculateEnergy(const Vec3 &position) const {
    auto difference = position - origin;
    scalar distanceFromOrigin = difference.norm();
    scalar distanceFromSphere = distanceFromOrigin - radius;
    auto energy = static_cast<scalar>(0.);
    if (distanceFromSphere < 0) {
        energy = static_cast<scalar>(0.5) * forceConstant * distanceFromSphere * distanceFromSphere;
    }
    return energy;
}

void SphereOut::calculateForce(Vec3 &force, const Vec3 &position) const {
    auto difference = position - origin;
    scalar distanceFromOrigin = difference.norm();
    scalar distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere < 0) {
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}

void SphereOut::calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const {
    auto difference = position - origin;
    scalar distanceFromOrigin = difference.norm();
    scalar distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere < 0) {
        energy += 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}

SphericalBarrier::SphericalBarrier(particle_type_type particleType, const Vec3 &origin, scalar radius, scalar height, scalar width)
        : super(particleType), origin(origin), radius(radius), height(height), width(width), r1(radius - width), r2(radius - width / static_cast<scalar>(2.)),
          r3(radius + width / static_cast<scalar>(2.)), r4(radius + width), effectiveForceConstant(static_cast<scalar>(4.) * height / width / width) {
    if (width > radius) {
        throw std::invalid_argument("SphericalBarrier must have a radius larger than its width");
    }
}

readdy::scalar SphericalBarrier::getRelevantLengthScale() const noexcept {
    return width;
}

readdy::scalar SphericalBarrier::getMaximalForce(readdy::scalar /*kbt*/) const noexcept {
    return static_cast<scalar>(2. * height / width);
}

readdy::scalar SphericalBarrier::calculateEnergy(const Vec3 &position) const {
    const auto difference = position - origin;
    const auto distance = difference.norm();
    if (distance < r1) {
        return static_cast<scalar>(0.);
    }
    if (r4 <= distance) {
        return static_cast<scalar>(0.);
    }
    if (r1 <= distance && distance < r2) {
        return static_cast<scalar>(0.5) * effectiveForceConstant * std::pow(distance - radius + width, static_cast<scalar>(2.));
    }
    if (r2 <= distance && distance < r3) {
        return height - static_cast<scalar>(0.5) * effectiveForceConstant * std::pow(distance - radius, static_cast<scalar>(2.));
    }
    if (r3 <= distance && distance < r4) {
        return static_cast<scalar>(0.5) * effectiveForceConstant * std::pow(distance - radius - width, static_cast<scalar>(2.));
    }
    throw std::runtime_error("Thou shalt not pass");
}

void SphericalBarrier::calculateForce(Vec3 &force, const Vec3 &position) const {
    const auto difference = position - origin;
    const auto distance = difference.norm();
    if (distance < r1) {
        // nothing happens
    } else if (r4 <= distance) {
        // nothing happens
    } else if (r1 <= distance && distance < r2) {
        force += - effectiveForceConstant * (distance - radius + width) * difference / distance;
    } else if (r2 <= distance && distance < r3) {
        force += effectiveForceConstant * (distance - radius) * difference / distance;
    } else if (r3 <= distance && distance < r4) {
        force += - effectiveForceConstant * (distance - radius - width) * difference / distance;
    } else {
        throw std::runtime_error("Not gonna happen");
    }
}

void SphericalBarrier::calculateForceAndEnergy(Vec3 &force, readdy::scalar &energy, const Vec3 &position) const {
    const auto difference = position - origin;
    const auto distance = difference.norm();
    if (distance < r1) {
        // nothing happens
    } else if (r4 <= distance) {
        // nothing happens
    } else if (r1 <= distance && distance < r2) {
        force += - effectiveForceConstant * (distance - radius + width) * difference / distance;
        energy += 0.5 * effectiveForceConstant * std::pow(distance - radius + width, 2.);
    } else if (r2 <= distance && distance < r3) {
        force += effectiveForceConstant * (distance - radius) * difference / distance;
        energy += height - 0.5 * effectiveForceConstant * std::pow(distance - radius, 2.);
    } else if (r3 <= distance && distance < r4) {
        force += - effectiveForceConstant * (distance - radius - width) * difference / distance;
        energy += 0.5 * effectiveForceConstant * std::pow(distance - radius - width, 2.);
    } else {
        throw std::runtime_error("Not gonna happen");
    }
}

std::string SphericalBarrier::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<SphericalBarrier>() << "[type: " << _particleType << ", origin: " << origin
       << ", radius: " << radius << ", height(energy): " << height << ", width: " << width << "]";
    return ss.str();
}

/////////////////////////////////////////////////////////////////////////////
//
// Potentials order 2
//
/////////////////////////////////////////////////////////////////////////////

/*
 * Harmonic repulsion
 */

HarmonicRepulsion::HarmonicRepulsion(particle_type_type type1, particle_type_type type2,
                                     scalar forceConstant, scalar interactionDistance)
        : super(type1, type2), _forceConstant(forceConstant), _interactionDistance(interactionDistance),
          _interactionDistanceSquared(interactionDistance*interactionDistance) {}

scalar HarmonicRepulsion::getForceConstant() const {
    return _forceConstant;
}

scalar HarmonicRepulsion::getMaximalForce(scalar /*kbt*/) const noexcept {
    return _forceConstant * getCutoffRadius();
}

scalar HarmonicRepulsion::calculateEnergy(const Vec3 &x_ij) const {
    auto distanceSquared = x_ij * x_ij;
    if (distanceSquared < _interactionDistanceSquared) {
        distanceSquared = std::sqrt(distanceSquared);
        distanceSquared -= _interactionDistance;
        distanceSquared *= distanceSquared;
        return static_cast<scalar>(0.5) * distanceSquared * getForceConstant();
    }
    return 0;
}

void HarmonicRepulsion::calculateForce(Vec3 &force, const Vec3 &x_ij) const {
    auto squared = x_ij * x_ij;
    if (squared < _interactionDistanceSquared && squared > 0) {
        squared = std::sqrt(squared);
        force = (getForceConstant() * (squared - _interactionDistance)) / squared * x_ij;
    } else {
        force = {0, 0, 0};
    }
}

void HarmonicRepulsion::calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &x_ij) const {
    auto squared = x_ij * x_ij;
    if (squared < _interactionDistanceSquared && squared > 0) {
        squared = std::sqrt(squared);
        energy += 0.5 * getForceConstant() * std::pow(squared - _interactionDistance, 2);
        force = (getForceConstant() * (squared - _interactionDistance)) / squared * x_ij;
    } else {
        force = {0, 0, 0};
    }
}

scalar HarmonicRepulsion::getCutoffRadius() const {
    return _interactionDistance;
}

scalar HarmonicRepulsion::getCutoffRadiusSquared() const {
    return _interactionDistanceSquared;
}

std::string HarmonicRepulsion::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<HarmonicRepulsion>() << "[type1: " << _particleType1 << ", type2: "
       << _particleType2 << ", forceConstant: " << _forceConstant << "]";
    return ss.str();
}

scalar HarmonicRepulsion::interactionDistance() const {
    return _interactionDistance;
}

/**
 * Weak interaction piecewise harmonic
 */

scalar WeakInteractionPiecewiseHarmonic::getMaximalForce(scalar /*kbt*/) const noexcept {
    scalar fMax1 = forceConstant * conf.desiredParticleDistance;
    scalar fMax2 = 2 * conf.depthAtDesiredDistance *
                   (conf.noInteractionDistance - conf.desiredParticleDistance);
    return std::max(fMax1, fMax2);
}


WeakInteractionPiecewiseHarmonic::WeakInteractionPiecewiseHarmonic(particle_type_type type1, particle_type_type type2,
                                                                   const scalar forceConstant,
                                                                   const Configuration &config)
        : super(type1, type2), forceConstant(forceConstant), conf(config) {}

std::ostream &operator<<(std::ostream &os, const WeakInteractionPiecewiseHarmonic::Configuration &configuration) {
    os << "desiredParticleDistance: " << configuration.desiredParticleDistance << " depthAtDesiredDistance: "
       << configuration.depthAtDesiredDistance << " noInteractionDistance: " << configuration.noInteractionDistance;
    return os;
}

scalar WeakInteractionPiecewiseHarmonic::calculateEnergy(const Vec3 &x_ij) const {
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

void WeakInteractionPiecewiseHarmonic::calculateForce(Vec3 &force, const Vec3 &x_ij) const {
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
        force = factor * x_ij / dist;
    } else {
        force = {0, 0, 0};
    }
}

void WeakInteractionPiecewiseHarmonic::calculateForceAndEnergy(Vec3 &force, scalar  &energy, const Vec3 &x_ij) const {
    energy += calculateEnergy(x_ij);
    calculateForce(force, x_ij);
}

scalar  WeakInteractionPiecewiseHarmonic::getCutoffRadius() const {
    return conf.noInteractionDistance;
}

scalar  WeakInteractionPiecewiseHarmonic::getCutoffRadiusSquared() const {
    return conf.noInteractionDistanceSquared;
}

std::string WeakInteractionPiecewiseHarmonic::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<HarmonicRepulsion>() << "[type1: " << _particleType1 << ", type2: "
       << _particleType2 << ", configuration[" << conf
       << "], forceConstant: " << forceConstant << "]";
    return ss.str();
}

WeakInteractionPiecewiseHarmonic::Configuration::Configuration(const scalar  desiredParticleDistance,
                                                               const scalar  depthAtDesiredDistance,
                                                               const scalar  noInteractionDistance)
        : desiredParticleDistance(desiredParticleDistance), depthAtDesiredDistance(depthAtDesiredDistance),
          noInteractionDistance(noInteractionDistance),
          noInteractionDistanceSquared(noInteractionDistance * noInteractionDistance) {}

LennardJones::LennardJones(particle_type_type type1, particle_type_type type2,
                           unsigned int m, unsigned int n, scalar  cutoffDistance,
                           bool shift, scalar  epsilon, scalar  sigma)
        : super(type1, type2), m(m), n(n),
          cutoffDistance(cutoffDistance), shift(shift), epsilon(epsilon), sigma(sigma),
          cutoffDistanceSquared(cutoffDistance * cutoffDistance) {
    if (m <= n) {
        throw std::invalid_argument("When constructing the LJ potential, the first exponent m=" + std::to_string(m) +
                                    " was not greater than the second exponent n=" + std::to_string(n) + "!");
    }
    auto dm = static_cast<scalar >(m);
    auto dn = static_cast<scalar >(n);
    scalar  r_min = sigma * std::pow(dn / dm, c_::one / (dn - dm));
    k = -epsilon / (std::pow(sigma / r_min, dm) - std::pow(sigma / r_min, dn));
}

scalar  LennardJones::calculateEnergy(const Vec3 &x_ij) const {
    const auto r = x_ij.norm();
    if (r > cutoffDistance) return 0;
    return shift ? energy(r) - energy(cutoffDistance) : energy(r);
}

void LennardJones::calculateForce(Vec3 &force, const Vec3 &x_ij) const {
    const auto norm = x_ij.norm();
    if(norm <= cutoffDistance) {
        force -= k * ( 1 / (sigma * sigma)) * (m * std::pow(sigma / norm, m + 2) - n * std::pow(sigma / norm, n + 2)) *
                 x_ij;
    }
}

void LennardJones::calculateForceAndEnergy(Vec3 &force, scalar  &energy, const Vec3 &x_ij) const {
    energy += calculateEnergy(x_ij);
    calculateForce(force, x_ij);
}

scalar LennardJones::getCutoffRadius() const {
    return cutoffDistance;
}

scalar LennardJones::getCutoffRadiusSquared() const {
    return cutoffDistanceSquared;
}

scalar LennardJones::energy(scalar  r) const {
    return k * (std::pow(sigma / r, m) - std::pow(sigma / r, n));
}

scalar LennardJones::getMaximalForce(scalar /*kbt*/) const noexcept {
    return 0;
}

std::string LennardJones::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<LennardJones>() << "[m: " << m << " n: "
       << n << " cutoffDistance: " << cutoffDistance << " shift: " << shift
       << " epsilon: " << epsilon << " k: " << k << "]";
    return ss.str();
}

LennardJones::~LennardJones() = default;


ScreenedElectrostatics::ScreenedElectrostatics(particle_type_type type1, particle_type_type type2,
                                               scalar  electrostaticStrength, scalar  inverseScreeningDepth,
                                               scalar  repulsionStrength, scalar  repulsionDistance, unsigned int exponent,
                                               scalar  cutoff)
        : super(type1, type2), electrostaticStrength(electrostaticStrength), inverseScreeningDepth(inverseScreeningDepth),
          repulsionStrength(repulsionStrength), repulsionDistance(repulsionDistance), exponent(exponent), cutoff(cutoff),
          cutoffSquared(cutoff * cutoff) {
    if (inverseScreeningDepth < 0) {
        throw std::invalid_argument("inverse screening depth must be positive!");
    }
    if (repulsionStrength < 0) {
        throw std::invalid_argument("repulsion strength must be positive!");
    }
    if (repulsionDistance < 0) {
        throw std::invalid_argument("repulsion distance must be positive!");
    }
    if (cutoff < 0) {
        throw std::invalid_argument("cutoff must be positive!");
    }
}

scalar ScreenedElectrostatics::getMaximalForce(scalar /*kbt*/) const noexcept {
    return 0;
}

scalar ScreenedElectrostatics::getCutoffRadius() const {
    return cutoff;
}

scalar ScreenedElectrostatics::getCutoffRadiusSquared() const {
    return cutoffSquared;
}

std::string ScreenedElectrostatics::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<ScreenedElectrostatics>() << "[electrostaticStrength: " << electrostaticStrength
       << " inverseScreeningDepth: " << inverseScreeningDepth
       << " repulsionStrength: " << repulsionStrength << " repulsionDistance: " << repulsionDistance
       << " exponent: " << exponent
       << " cutoff: " << cutoff << "]";
    return ss.str();
}

scalar  ScreenedElectrostatics::calculateEnergy(const Vec3 &x_ij) const {
    const scalar  distance = x_ij.norm();
    scalar  result = electrostaticStrength * std::exp(-inverseScreeningDepth * distance) / distance;
    result += repulsionStrength * std::pow(repulsionDistance / distance, exponent);
    return result;
}

void ScreenedElectrostatics::calculateForce(Vec3 &force, const Vec3 &x_ij) const {
    auto distance = x_ij.norm();
    auto forceFactor = electrostaticStrength * std::exp(-inverseScreeningDepth * distance);
    forceFactor *= (inverseScreeningDepth / distance + c_::one / std::pow(distance, c_::two));
    forceFactor += repulsionStrength * exponent / repulsionDistance * std::pow( repulsionDistance / distance, exponent + c_::one);
    force += forceFactor * (- c_::one * x_ij / distance);
}

void ScreenedElectrostatics::calculateForceAndEnergy(Vec3 &force, scalar  &energy, const Vec3 &x_ij) const {
    calculateForce(force, x_ij);
    energy += calculateEnergy(x_ij);
}

ScreenedElectrostatics::~ScreenedElectrostatics() = default;

}
}
}
