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


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <readdy/api/Simulation.h>
#include <readdy/io/File.h>
#include <readdy/common/nodelete.h>
#include "ExportSchemeApi.h"
#include "PyPotential.h"
#include "ExportObservables.h"

namespace py = pybind11;

using rvp = py::return_value_policy;
using sim = readdy::Simulation;
using kp = readdy::plugin::KernelProvider;
using vec = readdy::model::Vec3;
using pot2 = readdy::rpy::PotentialOrder2Wrapper;
using model = readdy::model::KernelStateModel;
using ctx = readdy::model::KernelContext;
using kern = readdy::model::Kernel;

void exportTopologies(py::module &);

// thin wrappers

void registerPotentialOrder2(sim &self, pot2 *potential) {
    self.registerPotentialOrder2(potential);
}

void setBoxSize(sim &self, const vec &size) { /* explicitly choose void(vec) signature */ self.setBoxSize(size); }

std::string getSelectedKernelType(sim &self) { /* discard const reference */ return self.getSelectedKernelType(); }

void addParticle(sim &self, const std::string &type, const vec &pos) { self.addParticle(pos[0], pos[1], pos[2], type); }


enum class ParticleTypeFlavor {
    NORMAL = 0, TOPOLOGY, MEMBRANE
};

// module
PYBIND11_PLUGIN (api) {

    py::module api("api", "ReaDDy c++-api python module");

    exportSchemeApi<readdy::api::ReaDDyScheme>(api, "ReaDDyScheme");
    exportSchemeApi<readdy::api::AdvancedScheme>(api, "AdvancedScheme");

    auto topologyModule = api.def_submodule("top");
    exportTopologies(topologyModule);

    py::enum_<ParticleTypeFlavor>(api, "ParticleTypeFlavor")
            .value("NORMAL", ParticleTypeFlavor::NORMAL)
            .value("TOPOLOGY", ParticleTypeFlavor::TOPOLOGY)
            .value("MEMBRANE", ParticleTypeFlavor::MEMBRANE);

    py::class_ <sim> simulation(api, "Simulation");
    simulation.def(py::init<>())
            .def_property("kbt", &sim::getKBT, &sim::setKBT)
            .def_property("periodic_boundary", &sim::getPeriodicBoundary, &sim::setPeriodicBoundary)
            .def_property("box_size", &sim::getBoxSize, &setBoxSize)
            .def("register_particle_type",
                 [](sim &self, const std::string &name, double diffusionCoefficient, double radius,
                    ParticleTypeFlavor flavor) {
                     readdy::model::Particle::flavor_t f = [=] {
                         switch (flavor) {
                             case ParticleTypeFlavor::NORMAL:
                                 return readdy::model::Particle::FLAVOR_NORMAL;
                             case ParticleTypeFlavor::MEMBRANE:
                                 return readdy::model::Particle::FLAVOR_MEMBRANE;
                             case ParticleTypeFlavor::TOPOLOGY:
                                 return readdy::model::Particle::FLAVOR_TOPOLOGY;
                         }
                     }();
                     return self.registerParticleType(name, diffusionCoefficient, radius);
                 }, py::arg("name"), py::arg("diffusion_coefficient"), py::arg("radius"),
                 py::arg("flavor") = ParticleTypeFlavor::NORMAL)
            .def("add_particle", [](sim &self, const std::string &type, const vec &pos) {
                self.addParticle(pos[0], pos[1], pos[2], type);
            })
            .def("is_kernel_selected", &sim::isKernelSelected)
            .def("get_selected_kernel_type", &getSelectedKernelType)
            .def("record_trajectory", &sim::recordTrajectory)
            .def("close_trajectory_file", &sim::closeTrajectoryFile)
            .def("register_potential_order_2", &registerPotentialOrder2)
            .def("register_potential_harmonic_repulsion", &sim::registerHarmonicRepulsionPotential)
            .def("register_potential_piecewise_weak_interaction",
                 &sim::registerWeakInteractionPiecewiseHarmonicPotential)
            .def("register_potential_box", &sim::registerBoxPotential)
            .def("register_potential_sphere_in", &sim::registerSphereInPotential)
            .def("register_potential_sphere_out", &sim::registerSphereOutPotential)
            .def("register_potential_lennard_jones", &sim::registerLennardJonesPotential)
            .def("register_potential_shielded_electrostatics", &sim::registerShieldedElectrostaticsPotential)
            .def("get_particle_positions", &sim::getParticlePositions)
            .def("register_reaction_conversion", &sim::registerConversionReaction, rvp::reference_internal)
            .def("register_reaction_enzymatic", &sim::registerEnzymaticReaction, rvp::reference_internal)
            .def("register_reaction_fission", &sim::registerFissionReaction, rvp::reference_internal)
            .def("register_reaction_fusion", &sim::registerFusionReaction, rvp::reference_internal)
            .def("register_reaction_decay", &sim::registerDecayReaction, rvp::reference_internal)
            .def("register_compartment_sphere", &sim::registerCompartmentSphere)
            .def("register_compartment_plane", &sim::registerCompartmentPlane)
            .def("get_recommended_time_step", &sim::getRecommendedTimeStep)
            .def("kernel_supports_topologies", &sim::kernelSupportsTopologies)
            .def("create_topology_particle", &sim::createTopologyParticle)
            .def("add_topology", &sim::addTopology, rvp::reference)
            .def("set_kernel", &sim::setKernel)
            .def("run_scheme_readdy", [](sim &self, bool defaults) {
                     return std::make_unique<readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>>(
                             self.runScheme<readdy::api::ReaDDyScheme>(defaults)
                     );
                 }, py::arg("defaults")
            )
            .def("run_scheme_advanced", [](sim &self, bool defaults) {
                     return std::make_unique<readdy::api::SchemeConfigurator<readdy::api::AdvancedScheme>>(
                             self.runScheme<readdy::api::AdvancedScheme>(defaults)
                     );
                 }, py::arg("defaults")
            )
            .def("run", [](sim &self, const readdy::time_step_type steps, const double timeStep) {
                py::gil_scoped_release release;
                self.run(steps, timeStep);
            });
    exportObservables(api, simulation);

    py::class_<kp, std::unique_ptr<kp, readdy::util::nodelete>>(api, "KernelProvider")
            .def_static("get", &kp::getInstance, rvp::reference)
            .def("load_from_dir", &kp::loadKernelsFromDirectory);

    py::class_<pot2>(api, "Pot2")
            .def(py::init<std::string, std::string, py::object, py::object>())
            .def("calc_energy", &pot2::calculateEnergy)
            .def("calc_force", &pot2::calculateForce);

    py::class_<kern>(api, "Kernel").def("get_name", &kern::getName, rvp::reference);

    return api.ptr();

}
