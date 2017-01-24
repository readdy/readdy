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

CubePotential::CubePotential(const std::string &particleType, double forceConstant, const Vec3 &origin,
                             const Vec3 &extent, bool considerParticleRadius)
        : super(particleType), origin(origin), extent(extent), forceConstant(forceConstant),
          considerParticleRadius(considerParticleRadius), min(getMinExtent(origin, extent)),
          max(getMaxExtent(origin, extent)), particleRadius(0) {}

const Vec3 &CubePotential::getOrigin() const { return origin; }

const Vec3 &CubePotential::getExtent() const { return extent; }

double CubePotential::getForceConstant() const { return forceConstant; }

bool CubePotential::isConsiderParticleRadius() const { return considerParticleRadius; }

double CubePotential::getParticleRadius() const { return particleRadius; }

double CubePotential::getMaximalForce(double kbt) const noexcept {
    return 0;
}

double CubePotential::getRelevantLengthScale() const noexcept {
    return std::min(extent[0], std::min(extent[1], extent[2]));
}

void CubePotential::configureForType(const KernelContext *const ctx, const unsigned int type) {
    particleRadius = ctx->getParticleRadius(type);
}

std::string CubePotential::describe() {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

std::ostream &operator<<(std::ostream &os, const CubePotential &potential) {
    os << getPotentialName<CubePotential>() << "[type: " << potential.particleType <<", origin: " << potential.origin << ", extent: "
       << potential.extent << ", min: " << potential.min << ", max: " << potential.max << ", forceConstant: "
       << potential.forceConstant << ", considerParticleRadius: " << potential.considerParticleRadius << "]";
    return os;
}

double CubePotential::calculateEnergy(const Vec3 &position) const {
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

void CubePotential::calculateForce(Vec3 &force, const Vec3 &position) const {
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

void CubePotential::calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const {
    energy += calculateEnergy(position);
    calculateForce(force, position);
}

/*
 * Sphere Potential
 */

const Vec3 &SpherePotential::getOrigin() const { return origin; }

double SpherePotential::getRadius() const { return radius; }

double SpherePotential::getForceConstant() const { return forceConstant; }

double SpherePotential::getRelevantLengthScale() const noexcept {
    return radius;
}

double SpherePotential::getMaximalForce(double kbt) const noexcept {
    return 0;
}

SpherePotential::SpherePotential(const std::string &particleType, double f, const Vec3 &origin, double radius)
        : super(particleType), origin(origin), radius(radius), forceConstant(f) {}

void SpherePotential::configureForType(const KernelContext *const ctx, const unsigned int type) {}

std::ostream &operator<<(std::ostream &os, const SpherePotential &potential) {
    os << getPotentialName<SpherePotential>() << "[type: " << potential.particleType << ", origin: " << potential.origin << ", radius: "
       << potential.radius << ", forceConstant: " << potential.forceConstant << "]";
    return os;
}

std::string SpherePotential::describe() {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

double SpherePotential::calculateEnergy(const Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    double energy = 0;
    if (distanceFromSphere > 0) {
        energy = 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
    }
    return energy;
}

void SpherePotential::calculateForce(Vec3 &force, const Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere > 0) {
        force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
    }
}

void SpherePotential::calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const {
    auto difference = position - origin;
    double distanceFromOrigin = sqrt(difference * difference);
    double distanceFromSphere = distanceFromOrigin - radius;
    if (distanceFromSphere > 0) {
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

double HarmonicRepulsion::getMaximalForce(double kbt) const noexcept {
    return forceConstant * getCutoffRadius();
}

void HarmonicRepulsion::configureForTypes(const KernelContext *const ctx, unsigned int type1, unsigned int type2) {
    auto r1 = ctx->getParticleRadius(type1);
    auto r2 = ctx->getParticleRadius(type2);
    sumOfParticleRadii = r1 + r2;
    sumOfParticleRadiiSquared = sumOfParticleRadii * sumOfParticleRadii;
}

std::ostream &operator<<(std::ostream &os, const HarmonicRepulsion &repulsion) {
    os << getPotentialName<HarmonicRepulsion>() << "[type1: " << repulsion.particleType1 <<", type2: " << repulsion.particleType2 << ", forceConstant: " << repulsion.forceConstant << "]";
    return os;
}

std::string HarmonicRepulsion::describe() {
    std::stringstream ss;
    ss << *this;
    return ss.str();
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


/**
 * Weak interaction piecewise harmonic
 */

double WeakInteractionPiecewiseHarmonic::getMaximalForce(double kbt) const noexcept {
    double fMax1 = forceConstant * conf.desiredParticleDistance;
    double fMax2 = 2 * conf.depthAtDesiredDistance *
                   (conf.noInteractionDistance - conf.desiredParticleDistance);
    return std::max(fMax1, fMax2);
}


void WeakInteractionPiecewiseHarmonic::configureForTypes(const KernelContext *const, unsigned int, unsigned int) {}

WeakInteractionPiecewiseHarmonic::WeakInteractionPiecewiseHarmonic(const std::string &particleType1,
                                                                   const std::string &particleType2,
                                                                   const double forceConstant,
                                                                   const Configuration &config)
        : super(particleType1, particleType2), forceConstant(forceConstant), conf(config) {}

std::ostream &operator<<(std::ostream &os, const WeakInteractionPiecewiseHarmonic &harmonic) {
    os << getPotentialName<HarmonicRepulsion>() << "[type1: " << harmonic.particleType1 <<", type2: " << harmonic.particleType2 << ", configuration[" << harmonic.conf
       << "], forceConstant: " << harmonic.forceConstant<< "]";
    return os;
}

std::ostream &operator<<(std::ostream &os, const WeakInteractionPiecewiseHarmonic::Configuration &configuration) {
    os << "desiredParticleDistance: " << configuration.desiredParticleDistance << " depthAtDesiredDistance: "
       << configuration.depthAtDesiredDistance << " noInteractionDistance: " << configuration.noInteractionDistance;
    return os;
}

std::string WeakInteractionPiecewiseHarmonic::describe() {
    std::stringstream ss;
    ss << *this;
    return ss.str();
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


WeakInteractionPiecewiseHarmonic::Configuration::Configuration(const double desiredParticleDistance,
                                                               const double depthAtDesiredDistance,
                                                               const double noInteractionDistance)
        : desiredParticleDistance(desiredParticleDistance), depthAtDesiredDistance(depthAtDesiredDistance),
          noInteractionDistance(noInteractionDistance),
          noInteractionDistanceSquared(noInteractionDistance * noInteractionDistance) {}
}
}
}

