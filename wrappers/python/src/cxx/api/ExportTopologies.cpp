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
 * @file ExportTopologies.cpp
 * @brief Impl file for exporting topology related classes and functionality
 * @author clonker
 * @date 04.02.17
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <readdy/model/Particle.h>
#include <readdy/model/topologies/GraphTopology.h>

namespace py = pybind11;
using rvp = py::return_value_policy;

using particle = readdy::model::Particle;
using topology_particle = readdy::model::TopologyParticle;
using topology = readdy::model::top::GraphTopology;
using topology_potential = readdy::model::top::TopologyPotential;
using bonded_potential = readdy::model::top::BondedPotential;
using angle_potential = readdy::model::top::AnglePotential;
using torsion_potential = readdy::model::top::TorsionPotential;
using harmonic_bond = readdy::model::top::HarmonicBondPotential;
using harmonic_angle = readdy::model::top::HarmonicAnglePotential;
using cosine_dihedral = readdy::model::top::CosineDihedralPotential;
using vec3 = readdy::model::Vec3;

void exportTopologies(py::module &m) {
    py::class_<topology_particle>(m, "TopologyParticle")
            .def("get_position", [](topology_particle& self) {return self.getPos();})
            .def("get_type", [](topology_particle& self) {return self.getType();})
            .def("get_id", [](topology_particle& self) {return self.getId();});
    py::class_<topology>(m, "Topology")
            .def("get_n_particles", &topology::getNParticles)
            .def("add_harmonic_angle_potential", [](topology &self, const harmonic_angle::angles_t &angles) {
                self.addAnglePotential<harmonic_angle>(angles);
            })
            .def("add_harmonic_bond_potential", [](topology &self, const harmonic_bond::bonds_t &bonds) {
                self.addBondedPotential<harmonic_bond>(bonds);
            })
            .def("add_cosine_dihedral_potential", [](topology &self, const cosine_dihedral::dihedrals_t &dihedrals) {
                self.addTorsionPotential<cosine_dihedral>(dihedrals);
            });
    py::class_<topology_potential>(m, "TopologyPotential");
    {
        py::class_<bonded_potential, topology_potential>(m, "BondedPotential");
        py::class_<harmonic_bond::bond_t>(m, "HarmonicBondPotentialBond")
                .def(py::init<std::size_t, std::size_t, double, double>())
                .def_readonly("idx1", &harmonic_bond::bond_t::idx1)
                .def_readonly("idx2", &harmonic_bond::bond_t::idx2)
                .def_readonly("length", &harmonic_bond::bond_t::length)
                .def_readonly("force_constant", &harmonic_bond::bond_t::forceConstant);
        py::class_<harmonic_bond, bonded_potential>(m, "HarmonicBondPotential")
                .def("get_bonds", &harmonic_bond::getBonds)
                .def("calculate_energy", &harmonic_bond::calculateEnergy)
                .def("calculate_force", [](harmonic_bond &self, const vec3 &x_ij, const harmonic_bond::bond_t &bond) {
                    vec3 force(0, 0, 0);
                    self.calculateForce(force, x_ij, bond);
                    return force;
                });
    }
    {
        py::class_<angle_potential, topology_potential>(m, "AnglePotential");
        py::class_<harmonic_angle::angle_t>(m, "HarmonicAnglePotentialAngle")
                .def(py::init<std::size_t, std::size_t, std::size_t, double, double>())
                .def_readonly("idx1", &harmonic_angle::angle_t::idx1)
                .def_readonly("idx2", &harmonic_angle::angle_t::idx2)
                .def_readonly("idx3", &harmonic_angle::angle_t::idx3)
                .def_readonly("equilibrium_angle", &harmonic_angle::angle_t::equilibriumAngle)
                .def_readonly("force_constant", &harmonic_angle::angle_t::forceConstant);
        py::class_<harmonic_angle, angle_potential>(m, "HarmonicAnglePotential")
                .def("get_angles", &harmonic_angle::getAngles)
                .def("calculate_energy", &harmonic_angle::calculateEnergy)
                .def("calculate_force", [](harmonic_angle &self, const vec3 &x_ij, const vec3 &x_kj,
                                           const harmonic_angle::angle_t &angle) {
                    vec3 f1(0, 0, 0), f2(0, 0, 0), f3(0, 0, 0);
                    self.calculateForce(f1, f2, f3, x_ij, x_kj, angle);
                    return std::make_tuple(std::move(f1), std::move(f2), std::move(f3));
                });
    }
    {
        py::class_<torsion_potential, topology_potential>(m, "TorsionPotential");
        py::class_<cosine_dihedral::dihedral_t>(m, "CosineDihedralPotentialDihedral")
                .def(py::init<std::size_t, std::size_t, std::size_t, std::size_t, double, double, double>())
                .def_readonly("idx1", &cosine_dihedral::dihedral_t::idx1)
                .def_readonly("idx2", &cosine_dihedral::dihedral_t::idx2)
                .def_readonly("idx3", &cosine_dihedral::dihedral_t::idx3)
                .def_readonly("idx4", &cosine_dihedral::dihedral_t::idx4)
                .def_readonly("force_constant", &cosine_dihedral::dihedral_t::forceConstant)
                .def_readonly("phi_0", &cosine_dihedral::dihedral_t::phi_0)
                .def_readonly("multiplicity", &cosine_dihedral::dihedral_t::multiplicity);
        py::class_<cosine_dihedral>(m, "CosineDihedralPotential")
                .def("get_dihedrals", &cosine_dihedral::getDihedrals)
                .def("calculate_energy", &cosine_dihedral::calculateEnergy)
                .def("calculate_force", [](cosine_dihedral &self, const vec3 &x_ji, const vec3 &x_kj, const vec3 &x_kl,
                                           const cosine_dihedral::dihedral_t &dih) {
                    vec3 f1(0, 0, 0), f2(0, 0, 0), f3(0, 0, 0), f4(0, 0, 0);
                    self.calculateForce(f1, f2, f3, f4, x_ji, x_kj, x_kl, dih);
                    return std::make_tuple(f1, f2, f3, f4);
                });
    }
}