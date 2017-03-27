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
 * Cube potential
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

Cube::Cube(const std::string &particleType, double forceConstant, const Vec3 &origin,
                             const Vec3 &extent, bool considerParticleRadius)
        : super(particleType), origin(origin), extent(extent), forceConstant(forceConstant),
          considerParticleRadius(considerParticleRadius), min(getMinExtent(origin, extent)),
          max(getMaxExtent(origin, extent)), particleRadius(0) {}

const Vec3 &Cube::getOrigin() const { return origin; }

const Vec3 &Cube::getExtent() const { return extent; }

double Cube::getForceConstant() const { return forceConstant; }

bool Cube::isConsiderParticleRadius() const { return considerParticleRadius; }

double Cube::getParticleRadius() const { return particleRadius; }

double Cube::getMaximalForce(double) const noexcept {
    return 0;
}

double Cube::getRelevantLengthScale() const noexcept {
    return std::min(extent[0], std::min(extent[1], extent[2]));
}

void Cube::configureForType(const KernelContext *const ctx, const particle_type_type type) {
    particleRadius = ctx->getParticleRadius(type);
}

double Cube::calculateEnergy(const Vec3 &position) const {
    auto r = particleRadius;
    if (!isConsiderParticleRadius()) r = 0;

    double energy = 0;

    for (auto i = 0; i < 3; ++i) {
        if (position[i] - r < min[i] || position[i] + r > max[i]) {
            if (position[i] - r < min[i]) {
                energy += 0.5 * forceConstant * (position[i] - r - min[i]) * (position[i] - r - min[i]);
            } else {
                energy += 0.5 * forceConstant * (position[i] + r - max[i]) * (position[i] + r - max[i]);
            }
        }
    }

    return energy;
}

void Cube::calculateForce(Vec3 &force, const Vec3 &position) const {
    auto r = particleRadius;
    if (!isConsiderParticleRadius()) r = 0;
    for (auto i = 0; i < 3; i++) {
        if (position[i] - r < min[i] || position[i] + r > max[i]) {
            if (position[i] - r < min[i]) {
                force[i] += -1 * forceConstant * (position[i] - r - min[i]);
            } else {
                force[i] += -1 * forceConstant * (position[i] + r - max[i]);
            }
        }
    }
}

void Cube::calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const {
    energy += calculateEnergy(position);
    calculateForce(force, position);
}

std::string Cube::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<Cube>() << "[type: " << particleType << ", origin: " << origin
       << ", extent: " << extent << ", min: " << min << ", max: " << max << ", forceConstant: "
       << forceConstant << ", considerParticleRadius: " << considerParticleRadius << "]";
    return ss.str();
}

/*
 * Sphere Potentials
 */

double SphereIn::getRelevantLengthScale() const noexcept {
    return radius;
}

double SphereIn::getMaximalForce(double) const noexcept {
    return 0;
}

SphereIn::SphereIn(const std::string &particleType, double f, const Vec3 &origin, double radius)
        : super(particleType), origin(origin), radius(radius), forceConstant(f) {}

void SphereIn::configureForType(const KernelContext *const, const particle_type_type) {}

double SphereIn::calculateEnergy(const Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    double energy = 0.;
    if (distanceFromSphere > 0) {
        energy = 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
    }
    return energy;
}

void SphereIn::calculateForce(Vec3 &force, const Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere > 0) {
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}

void SphereIn::calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere > 0) {
        energy += 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}

std::string SphereIn::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<SphereIn>() << "[type: " << particleType << ", origin: " << origin
       << ", radius: " << radius << ", forceConstant: " << forceConstant << "]";
    return ss.str();
}

SphereOut::SphereOut(const std::string &particleType, double forceConstant, const Vec3 &origin, double radius)
        : super(particleType), forceConstant(forceConstant), origin(origin), radius(radius) {}

void SphereOut::configureForType(const KernelContext *const ctx, const PotentialOrder1::particle_type_type type) {}

std::string SphereOut::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<SphereIn>() << "[type: " << particleType << ", origin: " << origin
       << ", radius: " << radius << ", forceConstant: " << forceConstant << "]";
    return ss.str();
}

double SphereOut::getRelevantLengthScale() const noexcept  {
    return radius;
}

double SphereOut::getMaximalForce(double kbt) const noexcept {
    return 0;
}

double SphereOut::calculateEnergy(const Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    double energy = 0.;
    if (distanceFromSphere < 0) {
        energy = 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
    }
    return energy;
}

void SphereOut::calculateForce(Vec3 &force, const Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere < 0) {
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}

void SphereOut::calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere < 0) {
        energy += 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}

/////////////////////////////////////////////////////////////////////////////
//
// Potentials order 2
//
/////////////////////////////////////////////////////////////////////////////

/*
 * Harmonic repulsion
 */

HarmonicRepulsion::HarmonicRepulsion(const std::string &type1, const std::string &type2, double forceConstant)
        : super(type1, type2), forceConstant(forceConstant), sumOfParticleRadii(-1), sumOfParticleRadiiSquared(-1) {}

double HarmonicRepulsion::getSumOfParticleRadii() const {
    return sumOfParticleRadii;
}

double HarmonicRepulsion::getSumOfParticleRadiiSquared() const {
    return sumOfParticleRadiiSquared;
}

double HarmonicRepulsion::getForceConstant() const {
    return forceConstant;
}

double HarmonicRepulsion::getMaximalForce(double) const noexcept {
    return forceConstant * getCutoffRadius();
}

void HarmonicRepulsion::configureForTypes(const KernelContext *const ctx, particle_type_type type1,
                                          particle_type_type type2) {
    auto r1 = ctx->getParticleRadius(type1);
    auto r2 = ctx->getParticleRadius(type2);
    sumOfParticleRadii = r1 + r2;
    sumOfParticleRadiiSquared = sumOfParticleRadii * sumOfParticleRadii;
}

double HarmonicRepulsion::calculateEnergy(const Vec3 &x_ij) const {
    auto distanceSquared = x_ij * x_ij;
    if (distanceSquared < getSumOfParticleRadiiSquared()) {
        distanceSquared = std::sqrt(distanceSquared);
        distanceSquared -= getSumOfParticleRadii();
        distanceSquared *= distanceSquared;
        return 0.5 * distanceSquared * getForceConstant();
    } else {
        return 0;
    }
}

void HarmonicRepulsion::calculateForce(Vec3 &force, const Vec3 &x_ij) const {
    auto squared = x_ij * x_ij;
    if (squared < getSumOfParticleRadiiSquared() && squared > 0) {
        squared = std::sqrt(squared);
        force = (getForceConstant() * (squared - getSumOfParticleRadii())) / squared * x_ij;
    } else {
        force = {0, 0, 0};
    }
}

void HarmonicRepulsion::calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &x_ij) const {
    auto squared = x_ij * x_ij;
    if (squared < getSumOfParticleRadiiSquared() && squared > 0) {
        squared = std::sqrt(squared);
        energy += 0.5 * getForceConstant() * std::pow(squared - getSumOfParticleRadii(), 2);
        force = (getForceConstant() * (squared - getSumOfParticleRadii())) / squared * x_ij;
    } else {
        force = {0, 0, 0};
    }
}

double HarmonicRepulsion::getCutoffRadius() const {
    return sumOfParticleRadii;
}

double HarmonicRepulsion::getCutoffRadiusSquared() const {
    return sumOfParticleRadiiSquared;
}

std::string HarmonicRepulsion::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<HarmonicRepulsion>() << "[type1: " << particleType1 << ", type2: "
       << particleType2 << ", forceConstant: " << forceConstant << "]";
    return ss.str();
}

/**
 * Weak interaction piecewise harmonic
 */

double WeakInteractionPiecewiseHarmonic::getMaximalForce(double) const noexcept {
    double fMax1 = forceConstant * conf.desiredParticleDistance;
    double fMax2 = 2 * conf.depthAtDesiredDistance *
                   (conf.noInteractionDistance - conf.desiredParticleDistance);
    return std::max(fMax1, fMax2);
}


void WeakInteractionPiecewiseHarmonic::configureForTypes(const KernelContext *const, particle_type_type,
                                                         particle_type_type) {}

WeakInteractionPiecewiseHarmonic::WeakInteractionPiecewiseHarmonic(const std::string &particleType1,
                                                                   const std::string &particleType2,
                                                                   const double forceConstant,
                                                                   const Configuration &config)
        : super(particleType1, particleType2), forceConstant(forceConstant), conf(config) {}

std::ostream &operator<<(std::ostream &os, const WeakInteractionPiecewiseHarmonic::Configuration &configuration) {
    os << "desiredParticleDistance: " << configuration.desiredParticleDistance << " depthAtDesiredDistance: "
       << configuration.depthAtDesiredDistance << " noInteractionDistance: " << configuration.noInteractionDistance;
    return os;
}

double WeakInteractionPiecewiseHarmonic::calculateEnergy(const Vec3 &x_ij) const {
    const auto dist = std::sqrt(x_ij * x_ij);
    const auto len_part2 = conf.noInteractionDistance - conf.desiredParticleDistance;
    if (dist < conf.desiredParticleDistance) {
        // repulsive as we are closer than the desired distance
        return .5 * forceConstant * (dist - conf.desiredParticleDistance) * (dist - conf.desiredParticleDistance) -
               conf.depthAtDesiredDistance;
    } else {
        // attractive as we are further (but not too far) apart than the desired distance
        if (dist < conf.desiredParticleDistance + .5 * len_part2) {
            return .5 * conf.depthAtDesiredDistance * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) *
                   (dist - conf.desiredParticleDistance) * (dist - conf.desiredParticleDistance) -
                   conf.depthAtDesiredDistance;
        } else {
            // if we are not too far apart but still further than in the previous case, attractive
            if (dist < conf.noInteractionDistance) {
                return -0.5 * conf.depthAtDesiredDistance * (1 / (0.5 * len_part2)) * (1 / (0.5 * len_part2)) *
                       (dist - conf.noInteractionDistance) * (dist - conf.noInteractionDistance);
            }
        }
    }
    return 0;
}

void WeakInteractionPiecewiseHarmonic::calculateForce(Vec3 &force, const Vec3 &x_ij) const {
    const auto dist = std::sqrt(x_ij * x_ij);
    const auto len_part2 = conf.noInteractionDistance - conf.desiredParticleDistance;
    double factor = 0;
    if (dist < conf.desiredParticleDistance) {
        // repulsive as we are closer than the desired distance
        factor = -1 * forceConstant * (conf.desiredParticleDistance - dist);
    } else {
        // attractive as we are further (but not too far) apart than the desired distance
        if (dist < conf.desiredParticleDistance + .5 * len_part2) {
            factor = -1 * conf.depthAtDesiredDistance * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) *
                     (conf.desiredParticleDistance - dist);
        } else {
            // if we are not too far apart but still further than in the previous case, attractive
            if (dist < conf.noInteractionDistance) {
                factor = conf.depthAtDesiredDistance * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) *
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

void WeakInteractionPiecewiseHarmonic::calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &x_ij) const {
    energy += calculateEnergy(x_ij);
    calculateForce(force, x_ij);
}

double WeakInteractionPiecewiseHarmonic::getCutoffRadius() const {
    return conf.noInteractionDistance;
}

double WeakInteractionPiecewiseHarmonic::getCutoffRadiusSquared() const {
    return conf.noInteractionDistanceSquared;
}

std::string WeakInteractionPiecewiseHarmonic::describe() const {
    std::ostringstream ss;
    ss << getPotentialName<HarmonicRepulsion>() << "[type1: " << particleType1 << ", type2: "
       << particleType2 << ", configuration[" << conf
       << "], forceConstant: " << forceConstant << "]";
    return ss.str();
}

WeakInteractionPiecewiseHarmonic::Configuration::Configuration(const double desiredParticleDistance,
                                                               const double depthAtDesiredDistance,
                                                               const double noInteractionDistance)
        : desiredParticleDistance(desiredParticleDistance), depthAtDesiredDistance(depthAtDesiredDistance),
          noInteractionDistance(noInteractionDistance),
          noInteractionDistanceSquared(noInteractionDistance * noInteractionDistance) {}

LennardJones::LennardJones(const std::string &particleType1, const std::string &particleType2,
                           unsigned int m, unsigned int n, double cutoffDistance,
                           bool shift, double epsilon, double sigma)
        : super(particleType1, particleType2), m(m), n(n),
          cutoffDistance(cutoffDistance), shift(shift), epsilon(epsilon), sigma(sigma),
          cutoffDistanceSquared(cutoffDistance * cutoffDistance) {
    if (m <= n) {
        throw std::invalid_argument("When constructing the LJ potential, the first exponent m=" + std::to_string(m) +
                                    " was not greater than the second exponent n=" + std::to_string(n) + "!");
    }
    double dm = static_cast<double>(m);
    double dn = static_cast<double>(n);
    double r_min = sigma * std::pow(dn / dm, 1. / (dn - dm));
    k = -epsilon / (std::pow(sigma / r_min, dm) - std::pow(sigma / r_min, dn));
}

void
LennardJones::configureForTypes(const KernelContext *const context, particle_type_type type1,
                                particle_type_type type2) {

}

double LennardJones::calculateEnergy(const Vec3 &x_ij) const {
    const auto r = x_ij.norm();
    if (r > cutoffDistance) return 0;
    else return shift ? energy(r) - energy(cutoffDistance) : energy(r);
}

void LennardJones::calculateForce(Vec3 &force, const Vec3 &x_ij) const {
    const auto norm = x_ij.norm();
    if(norm <= cutoffDistance) {
        force -= k * ( 1 / (sigma * sigma)) * (m * std::pow(sigma / norm, m + 2) - n * std::pow(sigma / norm, n + 2)) *
                 x_ij;
    }
}

void LennardJones::calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &x_ij) const {
    energy += calculateEnergy(x_ij);
    calculateForce(force, x_ij);
}

double LennardJones::getCutoffRadius() const {
    return cutoffDistance;
}

double LennardJones::getCutoffRadiusSquared() const {
    return cutoffDistanceSquared;
}

double LennardJones::energy(double r) const {
    return k * (std::pow(sigma / r, m) - std::pow(sigma / r, n));
}

double LennardJones::getMaximalForce(double kbt) const noexcept {
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


ScreenedElectrostatics::ScreenedElectrostatics(const std::string &particleType1, const std::string &particleType2,
                                               double electrostaticStrength, double inverseScreeningDepth,
                                               double repulsionStrength, double repulsionDistance, unsigned int exponent,
                                               double cutoff)
        : super(particleType1, particleType2), electrostaticStrength(electrostaticStrength), inverseScreeningDepth(inverseScreeningDepth),
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

void
ScreenedElectrostatics::configureForTypes(const KernelContext *const context, particle_type_type type1,
                                          particle_type_type type2) {

}

double ScreenedElectrostatics::getMaximalForce(double kbt) const noexcept {
    return 0;
}

double ScreenedElectrostatics::getCutoffRadius() const {
    return cutoff;
}

double ScreenedElectrostatics::getCutoffRadiusSquared() const {
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

double ScreenedElectrostatics::calculateEnergy(const Vec3 &x_ij) const {
    const double distance = x_ij.norm();
    double result = electrostaticStrength * std::exp(-inverseScreeningDepth * distance) / distance;
    result += repulsionStrength * std::pow(repulsionDistance / distance, exponent);
    return result;
}

void ScreenedElectrostatics::calculateForce(Vec3 &force, const Vec3 &x_ij) const {
    const double distance = x_ij.norm();
    double forceFactor = electrostaticStrength * std::exp(-inverseScreeningDepth * distance);
    forceFactor *= (inverseScreeningDepth / distance + 1. / std::pow(distance, 2));
    forceFactor += repulsionStrength * exponent / repulsionDistance * std::pow( repulsionDistance / distance, exponent + 1);
    force += forceFactor * (- 1. * x_ij / distance);
}

void ScreenedElectrostatics::calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &x_ij) const {
    calculateForce(force, x_ij);
    energy += calculateEnergy(x_ij);
}

ScreenedElectrostatics::~ScreenedElectrostatics() = default;

}
}
}
