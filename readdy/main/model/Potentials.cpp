/**
 * << detailed description >>
 *
 * @file Potentials.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */

#include <readdy/model/Kernel.h>
#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/potentials/PotentialsOrder2.h>

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

