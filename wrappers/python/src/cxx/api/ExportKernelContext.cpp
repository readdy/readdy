/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file ExportKernelContext.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.09.17
 * @copyright GPL-3
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <readdy/model/Context.h>
#include <readdy/common/boundary_condition_operations.h>

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
            .def("add_conversion", (ReactionRegistry::ReactionId(ReactionRegistry::*)(
    const std::string &, const std::string &, const std::string &, scalar)) &ReactionRegistry::addConversion)
    .def("add_enzymatic", (ReactionRegistry::ReactionId(ReactionRegistry::*)(
    const std::string &, const std::string &, const std::string &, const std::string &, scalar, scalar)) &ReactionRegistry::addEnzymatic)
    .def("add_fission", (ReactionRegistry::ReactionId(ReactionRegistry::*)(
    const std::string &, const std::string &, const std::string &, const std::string &, scalar, scalar, scalar, scalar)) &ReactionRegistry::addFission)
    .def("add_fusion", (ReactionRegistry::ReactionId(ReactionRegistry::*)(
    const std::string &, const std::string &, const std::string &, const std::string &, scalar, scalar, scalar, scalar)) &ReactionRegistry::addFusion)
    .def("add_decay", (ReactionRegistry::ReactionId(ReactionRegistry::*)(
    const std::string &, const std::string &, scalar)) &ReactionRegistry::addDecay);

    py::class_<ParticleTypeRegistry>(module, "ParticleTypeRegistry")
            .def("id_of", &ParticleTypeRegistry::idOf)
            .def("add", &ParticleTypeRegistry::add, "name"_a, "diffusion_constant"_a, "flavor"_a = 0)
            .def("diffusion_constant_of", [](const ParticleTypeRegistry &self, const std::string &type) {
                return self.diffusionConstantOf(type);
            })
            .def("set_diffusion_constant_of", [](ParticleTypeRegistry &self, const std::string &type, scalar value){
                self.diffusionConstantOf(type) = value;
            })
            .def("n_types", &ParticleTypeRegistry::nTypes)
            .def("name_of", &ParticleTypeRegistry::nameOf)
            .def_property_readonly("type_mapping", &ParticleTypeRegistry::typeMapping, rvp::reference_internal);

    py::class_<PotentialRegistry>(module, "PotentialRegistry")
            .def("add_box",
                 [](PotentialRegistry &self, const std::string &particleType, scalar forceConstant, const Vec3 &origin,
                    const Vec3 &extent) {
                     return self.addBox(particleType, forceConstant, origin, extent);
                 })
            .def("add_harmonic_repulsion",
                 [](PotentialRegistry &self, const std::string &type1, const std::string &type2, scalar forceConstant, scalar interactionDistance) {
                     return self.addHarmonicRepulsion(type1, type2, forceConstant, interactionDistance);
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
                 [](PotentialRegistry &self, const std::string &particleType, scalar height, scalar width,
                    const Vec3 &origin, scalar radius) {
                     return self.addSphericalBarrier(particleType, height, width, origin, radius);
                 })
            .def("add_external_order1", [](PotentialRegistry& self, readdy::model::potentials::PotentialOrder1* pot) {
                return self.addUserDefined(pot);
            }, py::keep_alive<1, 2>())
            .def("add_external_order2", [](PotentialRegistry& self, readdy::model::potentials::PotentialOrder2* pot) {
                return self.addUserDefined(pot);
            }, py::keep_alive<1, 2>());

    py::class_<readdy::model::potentials::PotentialOrder1>(module, "PotentialOrder1");
    py::class_<readdy::model::potentials::PotentialOrder2>(module, "PotentialOrder2");

    py::class_<readdy::api::Bond>(module, "BondedPotentialConfiguration")
            .def(py::init([](scalar forceConstant, scalar length, const std::string &type) {
                if(type != "harmonic") {
                    throw std::invalid_argument("only suppoted type: \"harmonic\"");
                }
                readdy::api::Bond bond {forceConstant, length, readdy::api::BondType::HARMONIC};
                return bond;
            }), "force_constant"_a, "length"_a, "type"_a="harmonic");

    py::class_<readdy::api::Angle>(module, "AnglePotentialConfiguration")
            .def(py::init([](scalar forceConstant, scalar equilibriumAngle, const std::string &type) {
                if(type != "harmonic") {
                    throw std::invalid_argument("only suppoted type: \"harmonic\"");
                }
                readdy::api::Angle angle {forceConstant, equilibriumAngle, readdy::api::AngleType::HARMONIC};
                return angle;
            }), "force_constant"_a, "equilibrium_angle"_a, "type"_a="harmonic");

    py::class_<readdy::api::TorsionAngle>(module, "TorsionPotentialConfiguration")
            .def(py::init([](scalar forceConstant, scalar multiplicity, scalar phi0, const std::string &type) {
                if(type != "cos_dihedral") {
                    throw std::invalid_argument("only supported type: \"cos_dihedral\"");
                }
                readdy::api::TorsionAngle angle {forceConstant, multiplicity, phi0, readdy::api::TorsionType::COS_DIHEDRAL};
                return angle;
            }), "force_constant"_a, "multiplicity"_a, "phi0"_a, "type"_a="cos_dihedral");

    py::class_<TopologyRegistry>(module, "TopologyRegistry")
            .def("add_type", [](TopologyRegistry &self, const std::string &type) { return self.addType(type); })
            .def("add_structural_reaction", [](TopologyRegistry &self, const std::string &type,
                                               const readdy::model::top::reactions::StructuralTopologyReaction &reaction) {
                self.addStructuralReaction(type, reaction);
            })
            .def("add_spatial_reaction",
                 [](TopologyRegistry &self, const std::string &descriptor, scalar rate, scalar radius) {
                     self.addSpatialReaction(descriptor, rate, radius);
                 })
            .def("configure_bond_potential", &TopologyRegistry::configureBondPotential)
            .def("configure_angle_potential", &TopologyRegistry::configureAnglePotential)
            .def("configure_torsion_potential", &TopologyRegistry::configureTorsionPotential);

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
            .def("describe", &KernelContext::describe)
            .def("validate", &KernelContext::validate)
            .def("bounding_box_vertices", &KernelContext::getBoxBoundingVertices)
            .def("calculate_max_cutoff", &KernelContext::calculateMaxCutoff)
            .def_property("record_reactions_with_positions",
                          [](const KernelContext &self) { return self.recordReactionsWithPositions(); },
                          [](KernelContext &self, bool value) { self.recordReactionsWithPositions() = value; })
            .def_property("record_reaction_counts",
                          [](const KernelContext &self) { return self.recordReactionCounts(); },
                          [](KernelContext &self, bool value) { self.recordReactionCounts() = value; })
            .def("set_kernel_configuration", &KernelContext::setKernelConfiguration)
            .def_property_readonly("particle_types", [](KernelContext &self) -> ParticleTypeRegistry&  { return self.particleTypes(); }, rvp::reference_internal)
            .def_property_readonly("reactions", [](KernelContext &self) -> ReactionRegistry& { return self.reactions(); }, rvp::reference_internal)
            .def_property_readonly("potentials", [](KernelContext &self) -> PotentialRegistry& { return self.potentials(); }, rvp::reference_internal)
            .def_property_readonly("topologies", [](KernelContext &self) -> TopologyRegistry& { return self.topologyRegistry(); }, rvp::reference_internal)
            .def_property_readonly("compartments", [](KernelContext &self) -> CompartmentRegistry&  { return self.compartments(); }, rvp::reference_internal)
            .def_property_readonly("shortest_difference_fun", [](const KernelContext &self) -> std::function<Vec3(const Vec3&, const Vec3 &)> {
                return [&self](const Vec3 &v1, const Vec3 &v2) {
                    return bcs::shortestDifference(v1, v2, self.boxSize(), self.periodicBoundaryConditions());
                };
            })
            .def_property_readonly("fix_position_fun", [](const KernelContext &self) -> std::function<void(Vec3 &)> {
                return [&self](Vec3 &position) {
                    return bcs::fixPosition(position, self.boxSize(), self.periodicBoundaryConditions());
                };
            })
            .def_property_readonly("dist_squared_fun", [](const KernelContext &self) -> std::function<scalar(const Vec3 &, const Vec3 &)> {
                return [&self](const Vec3 &p1, const Vec3 &p2) {
                    return bcs::distSquared(p1, p2, self.boxSize(), self.periodicBoundaryConditions());
                };
            });

}
