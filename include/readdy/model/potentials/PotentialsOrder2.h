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
