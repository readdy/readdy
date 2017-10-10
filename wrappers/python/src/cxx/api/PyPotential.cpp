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
 * @file PotentialWrapper.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.06.16
 */

#include "PyPotential.h"

namespace readdy {

namespace rpy {
PotentialOrder2Wrapper::PotentialOrder2Wrapper(particle_type_type particleType1, particle_type_type particleType2,
                                               pybind11::object o1, pybind11::object o2)
        : PotentialOrder2(particleType1, particleType2), calcEnergyFun(new pybind11::object(o1), [](pybind11::object *o) {
                      pybind11::gil_scoped_acquire lock;
                      delete o;
                  }), calcForceFun(new pybind11::object(o2), [](pybind11::object *o) {
                      pybind11::gil_scoped_acquire lock;
                      delete o;
                  }) {}

readdy::scalar PotentialOrder2Wrapper::calculateEnergy(const Vec3 &x_ij) const {
    pybind11::gil_scoped_acquire lock;
    return ((*calcEnergyFun)(x_ij)).cast<readdy::scalar>();
}

void PotentialOrder2Wrapper::calculateForce(Vec3 &force, const Vec3 &x_ij) const {
    pybind11::gil_scoped_acquire lock;
    force += ((*calcForceFun)(x_ij)).cast<Vec3>();
}

void
PotentialOrder2Wrapper::calculateForceAndEnergy(Vec3 &force, readdy::scalar &energy, const Vec3 &x_ij) const {
    pybind11::gil_scoped_acquire lock;
    energy += ((*calcEnergyFun)(x_ij)).cast<readdy::scalar>();
    force += ((*calcForceFun)(x_ij)).cast<Vec3>();
}

readdy::scalar PotentialOrder2Wrapper::getCutoffRadiusSquared() const {
    return getCutoffRadius() * getCutoffRadius();
}

std::string PotentialOrder2Wrapper::describe() const {
    return "Python wrapped potential order 2";
}

std::string PotentialOrder2Wrapper::type() const {
    return describe();
}


}


}

