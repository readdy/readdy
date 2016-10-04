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

#ifndef READDY_MAIN_POTENTIALSORDER2_H
#define READDY_MAIN_POTENTIALSORDER2_H

#include <readdy/model/KernelContext.h>
#include "PotentialOrder2.h"

namespace readdy {
namespace model {
class Kernel;
namespace potentials {

class HarmonicRepulsion : public PotentialOrder2 {

public:
    HarmonicRepulsion(const Kernel *const kernel);

    virtual HarmonicRepulsion *replicate() const override = 0;

    double getSumOfParticleRadii() const;

    double getSumOfParticleRadiiSquared() const;

    double getForceConstant() const;

    void setForceConstant(double forceConstant);

    virtual void configureForTypes(unsigned int type1, unsigned int type2) override;

    virtual double getMaximalForce(double kbt) const noexcept override;

protected:
    const Kernel *const kernel;
    double sumOfParticleRadii;
    double sumOfParticleRadiiSquared;
    double forceConstant = 0;

};

class WeakInteractionPiecewiseHarmonic : public PotentialOrder2 {

public:
    WeakInteractionPiecewiseHarmonic(const Kernel *const kernel);

    virtual WeakInteractionPiecewiseHarmonic *replicate() const override = 0;

    virtual void configureForTypes(unsigned int type1, unsigned int type2) override;


    void setDesiredParticleDistance(double desiredParticleDistance);

    void setForceConstant(double forceConstant);

    void setDepthAtDesiredDistance(double depthAtDesiredDistance);

    void setNoInteractionDistance(double noInteractionDistance);

    virtual double getMaximalForce(double kbt) const noexcept override;

protected:
    const Kernel *const kernel;
    double desiredParticleDistance;
    double forceConstant;
    double depthAtDesiredDistance;
    double noInteractionDistance, noInteractionDistanceSquared;
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
}
}
}

#endif //READDY_MAIN_POTENTIALSORDER2_H
