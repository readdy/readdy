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
 * SingleCPU declarations of potentials of second order (particle-particle interactions).
 *
 * @file PotentialsOrder2.h
 * @brief Declare potentials of second order.
 * @author clonker
 * @date 09.06.16
 */

#ifndef READDY_MAIN_POTENTIALSORDER2_H_H
#define READDY_MAIN_POTENTIALSORDER2_H_H

#include <readdy/model/potentials/PotentialsOrder2.h>
#include <readdy/model/Kernel.h>

namespace readdy {

namespace kernel {
namespace scpu {

namespace potentials {
class SCPUHarmonicRepulsion : public readdy::model::potentials::HarmonicRepulsion {
    using vec_t = readdy::model::Vec3;

public:
    SCPUHarmonicRepulsion(const readdy::model::Kernel *const kernel);

    virtual double calculateEnergy(const vec_t &x_ij) const override;

    virtual void calculateForce(vec_t &force, const vec_t &x_ij) const override;

    virtual void calculateForceAndEnergy(vec_t &force, double &energy, const vec_t &x_ij) const override;

    virtual double getCutoffRadius() const override;

    virtual double getCutoffRadiusSquared() const override;

};

class SCPUWeakInteractionPiecewiseHarmonic : public readdy::model::potentials::WeakInteractionPiecewiseHarmonic {
    using vec_t = readdy::model::Vec3;
public:
    SCPUWeakInteractionPiecewiseHarmonic(const readdy::model::Kernel *const kernel);

    virtual double calculateEnergy(const readdy::model::Vec3 &x_ij) const override;

    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &x_ij) const override;

    virtual void
    calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &x_ij) const override;

    virtual double getCutoffRadius() const override;

    virtual double getCutoffRadiusSquared() const override;

};

}
}
}
}
#endif //READDY_MAIN_POTENTIALSORDER2_H_H
