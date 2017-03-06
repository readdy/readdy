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
 * @file PotentialWrapper.h
 * @brief << brief description >>
 * @author clonker
 * @date 10.06.16
 */

#ifndef READDY_MAIN_POTENTIALWRAPPER_H
#define READDY_MAIN_POTENTIALWRAPPER_H

#include <pybind11/pybind11.h>
#include <readdy/model/potentials/PotentialOrder2.h>

namespace readdy {
namespace rpy {
class PotentialOrder2Wrapper : public readdy::model::potentials::PotentialOrder2 {
public:
    std::string describe() override;

    PotentialOrder2Wrapper(const std::string &particleType1, const std::string &particleType2,
                           pybind11::object o1, pybind11::object o2);

    virtual double calculateEnergy(const model::Vec3 &x_ij) const override;

    virtual void calculateForce(model::Vec3 &force, const model::Vec3 &x_ij) const override;

    virtual void calculateForceAndEnergy(model::Vec3 &force, double &energy, const model::Vec3 &x_ij) const override;

    virtual double getCutoffRadius() const override {
        // todo!
        return 50;
    }

    virtual double getMaximalForce(double kbt) const noexcept override {
        // todo!
        return 0;
    }

    virtual double getCutoffRadiusSquared() const override;


protected:
    void configureForTypes(const model::KernelContext *const context, particle_type_type type1, particle_type_type type2) override;

    std::shared_ptr<pybind11::object> calcEnergyFun;
    std::shared_ptr<pybind11::object> calcForceFun;
};
}
}

#endif //READDY_MAIN_POTENTIALWRAPPER_H
