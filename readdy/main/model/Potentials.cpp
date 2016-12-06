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

CubePotential::CubePotential(const Kernel *const kernel) : PotentialOrder1(getPotentialName<CubePotential>()),
                                                           kernel(kernel) {}

void CubePotential::configureForType(const unsigned int type) {
    particleRadius = kernel->getKernelContext().getParticleRadius(type);
    for (auto i = 0; i < 3; i++) {
        if (extent[i] > 0) {
            min[i] = origin[i];
            max[i] = origin[i] + extent[i];
        } else {
            min[i] = origin[i] + extent[i];
            max[i] = origin[i];
        }
    }
}

const Vec3 &CubePotential::getOrigin() const { return origin; }

void CubePotential::setOrigin(const Vec3 &origin) { CubePotential::origin = origin; }

const Vec3 &CubePotential::getExtent() const { return extent; }

void CubePotential::setExtent(const Vec3 &extent) { CubePotential::extent = extent; }

double CubePotential::getForceConstant() const { return forceConstant; }

void CubePotential::setForceConstant(double forceConstant) { CubePotential::forceConstant = forceConstant; }

bool CubePotential::isConsiderParticleRadius() const { return considerParticleRadius; }

void CubePotential::setConsiderParticleRadius(
        bool considerParticleRadius) { CubePotential::considerParticleRadius = considerParticleRadius; }

double CubePotential::getParticleRadius() const { return particleRadius; }

double CubePotential::getMaximalForce(double kbt) const noexcept {
    return 0;
}

double CubePotential::getRelevantLengthScale() const noexcept {
    return std::min(extent[0], std::min(extent[1], extent[2]));
}

/*
 * Sphere Potential
 */

SpherePotential::SpherePotential(const Kernel *kernel) : PotentialOrder1(getPotentialName<SpherePotential>()), kernel(kernel) {}

const Vec3 &SpherePotential::getOrigin() const { return origin; }

void SpherePotential::setOrigin(const Vec3 &origin) { SpherePotential::origin = origin; }

double SpherePotential::getRadius() const { return radius; }

void SpherePotential::setRadius(double radius) { SpherePotential::radius = radius; }

double SpherePotential::getForceConstant() const { return forceConstant; }

void SpherePotential::setForceConstant(double forceConstant) { SpherePotential::forceConstant = forceConstant; }

void SpherePotential::configureForType(const unsigned int type) {}

double SpherePotential::getRelevantLengthScale() const noexcept {
    return radius;
}

double SpherePotential::getMaximalForce(double kbt) const noexcept {
    return 0;
}

/////////////////////////////////////////////////////////////////////////////
//
// Potentials order 2
//
/////////////////////////////////////////////////////////////////////////////

/*
 * Harmonic repulsion
 */

void HarmonicRepulsion::configureForTypes(unsigned int type1, unsigned int type2) {
    auto r1 = kernel->getKernelContext().getParticleRadius(type1);
    auto r2 = kernel->getKernelContext().getParticleRadius(type2);
    sumOfParticleRadii = r1 + r2;
    sumOfParticleRadiiSquared = sumOfParticleRadii * sumOfParticleRadii;
}

HarmonicRepulsion::HarmonicRepulsion(const Kernel *const kernel) : PotentialOrder2(
        getPotentialName<HarmonicRepulsion>()), kernel(kernel) {}

double HarmonicRepulsion::getSumOfParticleRadii() const {
    return sumOfParticleRadii;
}

double HarmonicRepulsion::getSumOfParticleRadiiSquared() const {
    return sumOfParticleRadiiSquared;
}

double HarmonicRepulsion::getForceConstant() const {
    return forceConstant;
}

void HarmonicRepulsion::setForceConstant(double forceConstant) {
    HarmonicRepulsion::forceConstant = forceConstant;
}

double HarmonicRepulsion::getMaximalForce(double kbt) const noexcept {
    return forceConstant * getCutoffRadius();
}


/**
 * Weak interaction piecewise harmonic
 */

WeakInteractionPiecewiseHarmonic::WeakInteractionPiecewiseHarmonic(const readdy::model::Kernel *const kernel)
        : PotentialOrder2(getPotentialName<WeakInteractionPiecewiseHarmonic>()), kernel(kernel) {
    desiredParticleDistance = 0;
    forceConstant = 0;
    depthAtDesiredDistance = 0;
    noInteractionDistance = 0;
}

void WeakInteractionPiecewiseHarmonic::setDesiredParticleDistance(double desiredParticleDistance) {
    WeakInteractionPiecewiseHarmonic::desiredParticleDistance = desiredParticleDistance;
}

void WeakInteractionPiecewiseHarmonic::setForceConstant(double forceConstant) {
    WeakInteractionPiecewiseHarmonic::forceConstant = forceConstant;
}

void WeakInteractionPiecewiseHarmonic::setDepthAtDesiredDistance(double depthAtDesiredDistance) {
    WeakInteractionPiecewiseHarmonic::depthAtDesiredDistance = depthAtDesiredDistance;
}

void WeakInteractionPiecewiseHarmonic::setNoInteractionDistance(double noInteractionDistance) {
    WeakInteractionPiecewiseHarmonic::noInteractionDistance = noInteractionDistance;
    WeakInteractionPiecewiseHarmonic::noInteractionDistanceSquared = std::pow(noInteractionDistance, 2);
}

void WeakInteractionPiecewiseHarmonic::configureForTypes(unsigned int type1, unsigned int type2) {

}

double WeakInteractionPiecewiseHarmonic::getMaximalForce(double kbt) const noexcept {
    double fMax1 = forceConstant * desiredParticleDistance;
    double fMax2 = 2 * depthAtDesiredDistance * (noInteractionDistance - desiredParticleDistance);
    return std::max(fMax1, fMax2);
}


}
}
}

