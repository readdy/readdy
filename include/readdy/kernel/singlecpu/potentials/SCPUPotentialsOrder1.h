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
 * << detailed description >>
 *
 * @file P1Cube.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#ifndef READDY_MAIN_P1CUBE_H
#define READDY_MAIN_P1CUBE_H

#include <string>
#include <readdy/model/potentials/PotentialOrder1.h>
#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/Kernel.h>

namespace readdy {
namespace kernel {
namespace scpu {

namespace potentials {

class SCPUCubePotential : public readdy::model::potentials::CubePotential {
public:
    SCPUCubePotential(const readdy::model::Kernel *const kernel);

    virtual double calculateEnergy(const readdy::model::Vec3 &position) const override;

    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) const override;

    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy,
                                         const readdy::model::Vec3 &position) const override;

};

class SCPUSpherePotential : public readdy::model::potentials::SpherePotential {
public:
    SCPUSpherePotential(const readdy::model::Kernel * const kernel);

    virtual double calculateEnergy(const readdy::model::Vec3 &position) const override;

    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) const override;

    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &position) const override;

};

}
}
}

}

#endif //READDY_MAIN_P1CUBE_H
