/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file ExportKernelContext.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <readdy/model/Context.h>

namespace py = pybind11;

using rvp = py::return_value_policy;

using KernelContext = readdy::model::Context;
using ParticleTypeRegistry = readdy::model::ParticleTypeRegistry;
using ReactionRegistry = readdy::model::reactions::ReactionRegistry;
using PotentialRegistry = readdy::model::potentials::PotentialRegistry;
using TopologyRegistry = readdy::model::top::TopologyRegistry;
using CompartmentRegistry = readdy::model::compartments::CompartmentRegistry;

void exportKernelContext(py::module &module) {
    using namespace readdy;
    using namespace py::literals;

    py::class_<ReactionRegistry>(module, "ReactionRegistry")
            .def("add", &ReactionRegistry::add)
            .def("add_conversion", (ReactionRegistry::reaction_id(ReactionRegistry::*)(
    const std::string &, const std::string &, const std::string &, scalar)) &ReactionRegistry::addConversion)
    .def("add_enzymatic", (ReactionRegistry::reaction_id(ReactionRegistry::*)(
    const std::string &, const std::string &, const std::string &, const std::string &, scalar, scalar)) &ReactionRegistry::addEnzymatic)
    .def("add_fission", (ReactionRegistry::reaction_id(ReactionRegistry::*)(
    const std::string &, const std::string &, const std::string &, const std::string &, scalar, scalar, scalar, scalar)) &ReactionRegistry::addFission)
    .def("add_fusion", (ReactionRegistry::reaction_id(ReactionRegistry::*)(
    const std::string &, const std::string &, const std::string &, const std::string &, scalar, scalar, scalar, scalar)) &ReactionRegistry::addFusion)
    .def("add_decay", (ReactionRegistry::reaction_id(ReactionRegistry::*)(
    const std::string &, const std::string &, scalar)) &ReactionRegistry::addDecay);

    py::class_<ParticleTypeRegistry>(module, "ParticleTypeRegistry")
            .def("id_of", &ParticleTypeRegistry::id_of)
            .def("add", &ParticleTypeRegistry::add, "name"_a, "diffusion_constant"_a, "radius"_a, "flavor"_a = 0)
            .def("diffusion_constant_of", [](const ParticleTypeRegistry &self, const std::string &type) {
                return self.diffusion_constant_of(type);
            })
            .def("radius_of", [](const ParticleTypeRegistry &self, const std::string &type) {
                return self.radius_of(type);
            })
            .def("n_types", &ParticleTypeRegistry::n_types)
            .def("name_of", &ParticleTypeRegistry::name_of);

    py::class_<PotentialRegistry>(module, "PotentialRegistry")
            .def("add_box",
                 [](PotentialRegistry &self, const std::string &particleType, scalar forceConstant, const Vec3 &origin,
                    const Vec3 &extent, bool considerParticleRadius = true) {
                     return self.addCube(particleType, forceConstant, origin, extent, considerParticleRadius);
                 })
            .def("add_harmonic_repulsion",
                 [](PotentialRegistry &self, const std::string &type1, const std::string &type2, scalar forceConstant) {
                     return self.addHarmonicRepulsion(type1, type2, forceConstant);
                 })
            .def("add_weak_interaction_piecewise_harmonic",
                 [](PotentialRegistry &self, const std::string &type1, const std::string &type2,
                    scalar forceConstant, scalar desiredDist, scalar depth, scalar cutoff) {
                     return self.addWeakInteractionPiecewiseHarmonic(type1, type2, forceConstant, desiredDist, depth,
                                                                     cutoff);
                 })
            .def("add_lennard_jones",
                 [](PotentialRegistry &self, const std::string &type1, const std::string &type2, unsigned int m,
                    unsigned int n,
                    scalar cutoff, bool shift, scalar epsilon, scalar sigma) {
                     return self.addLennardJones(type1, type2, m, n, cutoff, shift, epsilon, sigma);
                 })
            .def("add_screened_electrostatics",
                 [](PotentialRegistry &self, const std::string &type1, const std::string &type2,
                    scalar electrostaticStrength, scalar inverseScreeningDepth,
                    scalar repulsionStrength, scalar repulsionDistance, unsigned int exponent,
                    scalar cutoff) {
                     return self.addScreenedElectrostatics(type1, type2, electrostaticStrength, inverseScreeningDepth,
                                                           repulsionStrength, repulsionDistance, exponent, cutoff);
                 })
            .def("add_sphere_out",
                 [](PotentialRegistry &self, const std::string &particleType, scalar forceConstant, const Vec3 &origin,
                    scalar radius) {
                     return self.addSphereOut(particleType, forceConstant, origin, radius);
                 })
            .def("add_sphere_in",
                 [](PotentialRegistry &self, const std::string &particleType, scalar forceConstant, const Vec3 &origin,
                    scalar radius) {
                     return self.addSphereIn(particleType, forceConstant, origin, radius);
                 })
            .def("add_spherical_barrier",
                 [](PotentialRegistry &self, const std::string &particleType, const Vec3 &origin, scalar radius,
                    scalar height,
                    scalar width) {
                     return self.addSphericalBarrier(particleType, origin, radius, height, width);
                 })
            .def("add_external_order1", [](PotentialRegistry& self, readdy::model::potentials::PotentialOrder1& pot) {
                return self.addUserDefined(&pot);
            })
            .def("add_external_order2", [](PotentialRegistry& self, readdy::model::potentials::PotentialOrder2& pot) {
                return self.addUserDefined(&pot);
            });

    py::class_<TopologyRegistry>(module, "TopologyRegistry")
            .def("add_type", [](TopologyRegistry &self, const std::string &type) { return self.add_type(type); })
            .def("add_structural_reaction", [](TopologyRegistry &self, const std::string &type,
                                               const readdy::model::top::reactions::StructuralTopologyReaction &reaction) {
                self.add_structural_reaction(type, reaction);
            })
            .def("add_spatial_reaction",
                 [](TopologyRegistry &self, const std::string &descriptor, scalar rate, scalar radius) {
                     self.add_spatial_reaction(descriptor, rate, radius);
                 })
            .def("configure_bond_potential", &TopologyRegistry::configure_bond_potential)
            .def("configure_angle_potential", &TopologyRegistry::configure_angle_potential)
            .def("configure_torsion_potential", &TopologyRegistry::configure_torsion_potential);

    py::class_<CompartmentRegistry>(module, "CompartmentRegistry")
            .def("add_sphere", [](CompartmentRegistry &self,
                                  const readdy::model::compartments::Compartment::label_conversion_map &conversions,
                                  const std::string &uniqueName,
                                  const Vec3 &origin, scalar radius, bool largerOrLess) {
                return self.addSphere(conversions, uniqueName, origin, radius, largerOrLess);
            })
            .def("add_plane", [](CompartmentRegistry &self,
                                 const readdy::model::compartments::Compartment::label_conversion_map &conversions,
                                 const std::string &uniqueName,
                                 const Vec3 &normalCoefficients, scalar distance, bool largerOrLess) {
                return self.addPlane(conversions, uniqueName, normalCoefficients, distance, largerOrLess);
            });

    py::class_<KernelContext>(module, "Context")
            .def(py::init<>())
            .def_property("kbt", [](const KernelContext &self) { return self.kBT(); },
                          [](KernelContext &self, scalar kbt) { self.kBT() = kbt; })
            .def("box_volume", &KernelContext::boxVolume)
            .def_property("box_size", [](const KernelContext &self) { return self.boxSize(); },
                          [](KernelContext &self, KernelContext::BoxSize boxSize) { self.boxSize() = boxSize; })
            .def_property("pbc", [](const KernelContext &self) { return self.periodicBoundaryConditions(); },
                          [](KernelContext &self, KernelContext::PeriodicBoundaryConditions pbc) {
                              self.periodicBoundaryConditions() = pbc;
                          })
            .def_property_readonly("dist_squared_fun", &KernelContext::distSquaredFun)
            .def_property_readonly("fix_position_fun", &KernelContext::fixPositionFun)
            .def_property_readonly("shortest_difference_fun", &KernelContext::shortestDifferenceFun)
            .def("configure", &KernelContext::configure, "debug_output"_a = true)
            .def("bounding_box_vertices", &KernelContext::getBoxBoundingVertices)
            .def("calculate_max_cutoff", &KernelContext::calculateMaxCutoff)
            .def_property("record_reactions_with_positions",
                          [](const KernelContext &self) { return self.recordReactionsWithPositions(); },
                          [](KernelContext &self, bool value) { self.recordReactionsWithPositions() = value; })
            .def_property("record_reaction_counts",
                          [](const KernelContext &self) { return self.recordReactionCounts(); },
                          [](KernelContext &self, bool value) { self.recordReactionCounts() = value; })
            .def("set_kernel_configuration", &KernelContext::setKernelConfiguration)
            .def("particle_types", [](KernelContext &self) -> ParticleTypeRegistry&  { return self.particle_types(); }, rvp::reference_internal)
            .def("reactions", [](KernelContext &self) -> ReactionRegistry& { return self.reactions(); }, rvp::reference_internal)
            .def("potentials", [](KernelContext &self) -> PotentialRegistry& { return self.potentials(); }, rvp::reference_internal)
            .def("topologies", [](KernelContext &self) -> TopologyRegistry& { return self.topology_registry(); }, rvp::reference_internal)
            .def("compartments", [](KernelContext &self) -> CompartmentRegistry&  { return self.compartments(); }, rvp::reference_internal);

}