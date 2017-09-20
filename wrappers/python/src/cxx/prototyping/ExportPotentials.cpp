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


//
// Created by mho on 10/08/16.
//
#include <pybind11/pybind11.h>

#include <readdy/model/potentials/PotentialOrder1.h>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/potentials/PotentialsOrder2.h>

namespace py = pybind11;
namespace pot = readdy::model::potentials;

using rvp = py::return_value_policy;

using rdy_pot = pot::Potential;
using rdy_pot1 = pot::PotentialOrder1;
using rdy_pot2 = pot::PotentialOrder2;

class PyPotentialO1 : public rdy_pot1 {
public:
    using rdy_pot1::PotentialOrder1;

    std::string describe() const override {
        return "User defined potential order 1 for type " + particleType;
    }

    virtual readdy::scalar calculateEnergy(const readdy::Vec3 &position) const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(readdy::scalar, rdy_pot1, "calculate_energy", calculateEnergy, position);
    }

    virtual readdy::Vec3 calculateForceInternal(const readdy::Vec3 &pos) const {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(readdy::Vec3, PyPotentialO1, "calculate_force", calculateForce, pos);
    }

    virtual void calculateForce(readdy::Vec3 &force, const readdy::Vec3 &position) const override {
        force += calculateForceInternal(position);
    }

    virtual void calculateForceAndEnergy(readdy::Vec3 &force, readdy::scalar &energy, const readdy::Vec3 &position) const override {
        calculateForce(force, position);
        energy += calculateEnergy(position);
    }

    virtual readdy::scalar getRelevantLengthScale() const noexcept override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(readdy::scalar, rdy_pot1, "get_relevant_length_scale", getRelevantLengthScale,);
    }

    virtual readdy::scalar getMaximalForce(readdy::scalar kbt) const noexcept override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(readdy::scalar, rdy_pot1, "get_maximal_force", getMaximalForce, kbt);
    }

protected:
    void configureForType(const readdy::model::ParticleTypeRegistry *const context, const particle_type_type type) override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(void, rdy_pot1, "configure_for_type", configureForType, type);
    }

};

class PyPotentialO2 : public rdy_pot2 {
public:

    using rdy_pot2::PotentialOrder2;

    std::string describe() const override {
        return "User defined potential for types " + std::to_string(particleType1) + " and " + std::to_string(particleType2);
    }

    virtual readdy::scalar getMaximalForce(readdy::scalar kbt) const noexcept override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(readdy::scalar, rdy_pot2, "get_maximal_force", getMaximalForce, kbt);
    }

    virtual readdy::scalar calculateEnergy(const readdy::Vec3 &x_ij) const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(readdy::scalar, rdy_pot2, "calculate_energy", calculateEnergy, x_ij);
    }

    virtual readdy::Vec3 calculateForceInternal(const readdy::Vec3 &pos) const {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(readdy::Vec3, PyPotentialO2, "calculate_force", calculateForce, pos);
    }

    virtual void calculateForce(readdy::Vec3 &force, const readdy::Vec3 &x_ij) const override {
        force += calculateForceInternal(x_ij);
    }

    virtual void calculateForceAndEnergy(readdy::Vec3 &force, readdy::scalar &energy, const readdy::Vec3 &x_ij) const override {
        calculateForce(force, x_ij);
        energy += calculateEnergy(x_ij);
    }

    virtual readdy::scalar getCutoffRadiusSquared() const override {
        const auto cutoff = getCutoffRadius();
        return cutoff * cutoff;
    }

    virtual readdy::scalar getCutoffRadius() const override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(readdy::scalar, rdy_pot2, "get_cutoff_radius", getCutoffRadius);
    }

protected:
    void configureForTypes(const readdy::model::ParticleTypeRegistry *const context, particle_type_type type1,
                           particle_type_type type2) override {
        py::gil_scoped_acquire gil;
        PYBIND11_OVERLOAD_PURE_NAME(void, rdy_pot2, "configure_for_types", configureForType, type1, type2);
    }


};

void exportPotentials(py::module &proto) {

    py::class_<rdy_pot>(proto, "Potential");
    py::class_<pot::Cube, rdy_pot>(proto, pot::getPotentialName<pot::Cube>().c_str())
            .def("get_name", &pot::getPotentialName < pot::Cube> );
    py::class_<pot::HarmonicRepulsion, rdy_pot>(proto, pot::getPotentialName<pot::HarmonicRepulsion>().c_str())
            .def("get_name", &pot::getPotentialName < pot::HarmonicRepulsion > );
    py::class_<pot::WeakInteractionPiecewiseHarmonic, rdy_pot>(proto,
                                                               pot::getPotentialName<pot::WeakInteractionPiecewiseHarmonic>().c_str())
            .def("get_name", &pot::getPotentialName < pot::WeakInteractionPiecewiseHarmonic > );

    py::class_<rdy_pot1, PyPotentialO1>(proto, "PotentialOrder1")
            .def(py::init<readdy::particle_type_type>())
            .def("calculate_energy", &rdy_pot1::calculateEnergy)
            .def("calculate_force", &rdy_pot1::calculateForce)
            .def("get_relevant_length_scale", &rdy_pot1::getRelevantLengthScale)
            .def("get_maximal_force", &rdy_pot1::getMaximalForce);
    py::class_<rdy_pot2, PyPotentialO2>(proto, "PotentialOrder2")
            .def(py::init<readdy::particle_type_type, readdy::particle_type_type>())
            .def("calculate_energy", &rdy_pot2::calculateEnergy)
            .def("calculate_force", &rdy_pot2::calculateForce)
            .def("get_cutoff_radius", &rdy_pot2::getCutoffRadius)
            .def("get_maximal_force", &rdy_pot2::getMaximalForce);
}