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

#include <json.hpp>

#include <readdy/api/Simulation.h>
#include <readdy/common/nodelete.h>
#include "ExportSchemeApi.h"
#include "PyPotential.h"
#include "ExportObservables.h"

namespace py = pybind11;

using rvp = py::return_value_policy;
using sim = readdy::Simulation;
using kp = readdy::plugin::KernelProvider;
using vec = readdy::Vec3;
using pot2 = readdy::rpy::PotentialOrder2Wrapper;
using model = readdy::model::StateModel;
using ctx = readdy::model::Context;
using kern = readdy::model::Kernel;

void exportTopologies(py::module &);

void exportKernelContext(py::module &);

void registerPotentialOrder2(sim &self, pot2 *potential) {
    self.registerPotentialOrder2(potential);
}

// thin wrappers
void setBoxSize(sim &self, const vec &size) { /* explicitly choose void(vec) signature */ self.setBoxSize(size); }

std::string getSelectedKernelType(sim &self) { /* discard const reference */ return self.getSelectedKernelType(); }

void addParticle(sim &self, const std::string &type, const vec &pos) { self.addParticle(type, pos[0], pos[1], pos[2]); }


enum class ParticleTypeFlavor {
    NORMAL = 0, TOPOLOGY = 1, MEMBRANE = 2
};

void exportApi(py::module &api) {
    using namespace pybind11::literals;
    exportSchemeApi<readdy::api::ReaDDyScheme>(api, "ReaDDyScheme");
    exportSchemeApi<readdy::api::AdvancedScheme>(api, "AdvancedScheme");

    auto topologyModule = api.def_submodule("top");
    exportTopologies(topologyModule);

    exportKernelContext(api);

    py::enum_<ParticleTypeFlavor>(api, "ParticleTypeFlavor")
            .value("NORMAL", ParticleTypeFlavor::NORMAL)
            .value("TOPOLOGY", ParticleTypeFlavor::TOPOLOGY)
            .value("MEMBRANE", ParticleTypeFlavor::MEMBRANE);

    py::enum_<readdy::api::BondType>(api, "BondType").value("HARMONIC", readdy::api::BondType::HARMONIC);
    py::enum_<readdy::api::AngleType>(api, "AngleType").value("HARMONIC", readdy::api::AngleType::HARMONIC);
    py::enum_<readdy::api::TorsionType>(api, "TorsionType").value("COS_DIHEDRAL", readdy::api::TorsionType::COS_DIHEDRAL);

    py::class_<sim> simulation(api, "Simulation");
    simulation.def(py::init<>())
            .def_property("kbt", &sim::getKBT, &sim::setKBT)
            .def_property("periodic_boundary", &sim::getPeriodicBoundary, &sim::setPeriodicBoundary)
            .def_property("box_size", &sim::getBoxSize, &setBoxSize)
            .def_property_readonly("single_precision", &sim::singlePrecision)
            .def_property_readonly("double_precision", &sim::doublePrecision)
            .def("register_particle_type",
                 [](sim &self, const std::string &name, readdy::scalar diffusionCoefficient,
                    ParticleTypeFlavor flavor) {
                     readdy::model::particle_flavor f = [=] {
                         switch (flavor) {
                             case ParticleTypeFlavor::NORMAL:
                                 return readdy::model::particleflavor::NORMAL;
                             case ParticleTypeFlavor::MEMBRANE:
                                 return readdy::model::particleflavor::MEMBRANE;
                             case ParticleTypeFlavor::TOPOLOGY:
                                 return readdy::model::particleflavor::TOPOLOGY;
                         }
                     }();
                     return self.registerParticleType(name, diffusionCoefficient, f);
                 }, "name"_a, "diffusion_coefficient"_a, "flavor"_a = ParticleTypeFlavor::NORMAL)
            .def("add_particle", [](sim &self, const std::string &type, const vec &pos) {
                self.addParticle(type, pos[0], pos[1], pos[2]);
            }, "type"_a, "pos"_a)
            .def("add_particles", [](sim &self, const std::string &type, const py::array_t<readdy::scalar> &particles) {
                auto nParticles = particles.shape(0);
                for(std::size_t i = 0; i < nParticles; ++i) {
                    self.addParticle(type, particles.at(i, 0), particles.at(i, 1), particles.at(i, 2));
                }
            })
            .def("set_kernel_config", &sim::setKernelConfiguration)
            .def("is_kernel_selected", &sim::isKernelSelected)
            .def("get_selected_kernel_type", &getSelectedKernelType)
            .def("register_potential_order_2", &registerPotentialOrder2, "potential"_a)
            .def("register_potential_harmonic_repulsion", &sim::registerHarmonicRepulsionPotential,
                 "type_a"_a, "type_b"_a, "force_constant"_a, "interaction_distance"_a)
            .def("register_potential_piecewise_weak_interaction",
                 &sim::registerWeakInteractionPiecewiseHarmonicPotential, "type_a"_a, "type_b"_a, "force_constant"_a,
                 "desired_particle_distance"_a, "depth"_a, "no_interaction_distance"_a)
            .def("register_potential_box", &sim::registerBoxPotential, "particle_type"_a, "force_constant"_a,
                 "origin"_a, "extent"_a)
            .def("register_potential_sphere_in", &sim::registerSphereInPotential, "particle_type"_a, "force_constant"_a,
                 "origin"_a, "radius"_a)
            .def("register_potential_sphere_out", &sim::registerSphereOutPotential, "particle_type"_a,
                 "force_constant"_a, "origin"_a, "radius"_a)
            .def("register_potential_spherical_barrier", &sim::registerSphericalBarrier, "particle_type"_a, "height"_a,
                 "width"_a, "origin"_a, "radius"_a)
            .def("register_potential_lennard_jones", &sim::registerLennardJonesPotential, "particle_type_a"_a,
                 "particle_type_b"_a, "exponent_m"_a, "exponent_n"_a, "cutoff"_a, "shift"_a, "epsilon"_a, "sigma"_a)
            .def("register_potential_screened_electrostatics", &sim::registerScreenedElectrostaticsPotential,
                 "particle_type_a"_a, "particle_type_b"_a, "electrostatic_strength"_a, "inverse_screening_depth"_a,
                 "repulsion_strength"_a, "repulsion_distance"_a, "exponent"_a, "cutoff"_a)
            .def("get_particle_positions", &sim::getParticlePositions)
            .def("register_reaction_conversion", &sim::registerConversionReaction, rvp::reference_internal,
                 "label"_a, "from_type"_a, "to_type"_a, "rate"_a)
            .def("register_reaction_enzymatic", &sim::registerEnzymaticReaction, rvp::reference_internal,
                 "label"_a, "catalyst_type"_a, "from_type"_a, "to_type"_a, "rate"_a, "educt_distance"_a)
            .def("register_reaction_fission", &sim::registerFissionReaction, rvp::reference_internal,
                 "label"_a, "from_type"_a, "to_type1"_a, "to_type2"_a, "rate"_a, "product_distance"_a,
                 "weight1"_a = .5, "weight2"_a = .5)
            .def("register_reaction_fusion", &sim::registerFusionReaction, rvp::reference_internal, "label"_a,
                 "from_type1"_a, "from_type2"_a, "to_type"_a, "rate"_a, "educt_distance"_a,
                 "weight1"_a = .5, "weight2"_a = .5)
            .def("register_reaction_decay", &sim::registerDecayReaction, rvp::reference_internal,
                 "label"_a, "particle_type"_a, "rate"_a)
            .def("register_spatial_topology_reaction", &sim::registerSpatialTopologyReaction,
                 "descriptor"_a, "rate"_a, "radius"_a)
            .def("register_compartment_sphere", &sim::registerCompartmentSphere,
                 "conversion_map"_a, "name"_a, "origin"_a, "radius"_a, "larger_or_less"_a)
            .def("register_compartment_plane", &sim::registerCompartmentPlane, "conversion_map"_a, "name"_a,
                 "normal_coefficients"_a, "distance_from_plane"_a, "larger_or_less"_a)
            .def("kernel_supports_topologies", &sim::kernelSupportsTopologies)
            .def("create_topology_particle", &sim::createTopologyParticle, "type"_a, "position"_a)
            .def("configure_topology_bond_potential", &sim::configureTopologyBondPotential, "type1"_a, "type2"_a,
                 "force_constant"_a, "length"_a, "type"_a = readdy::api::BondType::HARMONIC)
            .def("configure_topology_angle_potential", &sim::configureTopologyAnglePotential, "type1"_a, "type2"_a,
                 "type3"_a, "force_constant"_a, "equilibrium_angle"_a, "type"_a = readdy::api::AngleType::HARMONIC)
            .def("configure_topology_dihedral_potential", &sim::configureTopologyTorsionPotential, "type1"_a,
                 "type2"_a, "type3"_a, "type4"_a, "force_constant"_a, "multiplicity"_a, "phi_0"_a,
                 "type"_a = readdy::api::TorsionType::COS_DIHEDRAL)
            .def("register_topology_type", [](sim &self, const std::string& name) {
                return self.registerTopologyType(name);
            })
            .def("register_structural_topology_reaction", &sim::registerStructuralTopologyReaction)
            .def("get_particles_for_topology", &sim::getParticlesForTopology, "topology"_a)
            .def("add_topology", [](sim &self, const std::string &name,
                                    const std::vector<readdy::model::TopologyParticle> &particles) {
                return self.addTopology(name, particles);
            }, rvp::reference_internal, "type"_a, "particles"_a)
            .def("add_topology", [](sim &self, const std::string &name, const std::vector<std::string> &types,
                                    const py::array_t<readdy::scalar> &positions) {
                auto nParticles = positions.shape(0);
                auto nTypes = types.size();
                if(nParticles != nTypes && nTypes != 1) {
                    throw std::invalid_argument(fmt::format("the number of particles ({}) must be equal to the "
                                                                    "number of types ({})!", nParticles, nTypes));
                }
                std::vector<readdy::model::TopologyParticle> particles;
                for(std::size_t i = 0; i < nParticles; ++i) {
                    auto type =  nTypes != 1 ? types[i] : types[0];
                    particles.push_back(self.createTopologyParticle(type, readdy::Vec3(positions.at(i, 0),
                                                                                       positions.at(i, 1),
                                                                                       positions.at(i, 2))));
                }
                return self.addTopology(name, particles);
            }, rvp::reference_internal)
            .def("current_topologies", &sim::currentTopologies)
            .def("set_kernel", static_cast<void (sim::*)(const std::string&)>(&sim::setKernel), "name"_a)
            .def_property("context", [](sim &self) -> readdy::model::Context & {
                return self.currentContext();
            }, [](sim &self, const readdy::model::Context &context) {
                self.currentContext() = context;
            })
            .def("run_scheme_readdy", [](sim &self, bool defaults) {
                     return std::make_unique<readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>>(
                             self.runScheme<readdy::api::ReaDDyScheme>(defaults)
                     );
                 }, "defaults"_a = true
            )
            .def("run_scheme_advanced", [](sim &self, bool defaults) {
                     return std::make_unique<readdy::api::SchemeConfigurator<readdy::api::AdvancedScheme>>(
                             self.runScheme<readdy::api::AdvancedScheme>(defaults)
                     );
                 }, "defaults"_a = true
            )
            .def("run", [](sim &self, const readdy::time_step_type steps, const readdy::scalar timeStep) {
                py::gil_scoped_release release;
                self.run(steps, timeStep);
            }, "n_steps"_a, "time_step"_a)
            .def("performance_root", &sim::performanceRoot, rvp::reference);
    exportObservables(api, simulation);

    py::class_<kp, std::unique_ptr<kp, readdy::util::nodelete>>(api, "KernelProvider")
            .def_static("get", &kp::getInstance, rvp::reference)
            .def("load_from_dir", &kp::loadKernelsFromDirectory, "directory"_a)
            .def("available_kernels", &kp::availableKernels);

    py::class_<pot2>(api, "Pot2")
            .def(py::init<readdy::particle_type_type, readdy::particle_type_type, py::object, py::object>())
            .def("calc_energy", &pot2::calculateEnergy, "x_ij"_a)
            .def("calc_force", &pot2::calculateForce, "force"_a, "x_ij"_a);

    py::class_<kern>(api, "Kernel").def("get_name", &kern::getName, rvp::reference);
}
