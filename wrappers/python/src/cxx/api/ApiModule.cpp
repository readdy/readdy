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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <json.hpp>

#include <readdy/api/Simulation.h>
#include <readdy/model/actions/UserDefinedAction.h>
#include "ExportObservables.h"
#include "../common/ReadableParticle.h"
#include "PyTopology.h"

namespace py = pybind11;
namespace rma = readdy::model::actions;

using rvp = py::return_value_policy;
using sim = readdy::Simulation;
using kp = readdy::plugin::KernelProvider;
using vec = readdy::Vec3;
using model = readdy::model::StateModel;
using ctx = readdy::model::Context;
using kern = readdy::model::Kernel;

void exportTopologies(py::module &);
void exportLoopApi(py::module &);
void exportKernelContext(py::module &);

std::string getSelectedKernelType(sim &self) { /* discard const reference */ return self.selectedKernelType(); }

void addParticle(sim &self, const std::string &type, const vec &pos) { self.addParticle(type, pos[0], pos[1], pos[2]); }


enum class ParticleTypeFlavor {
    NORMAL = 0, TOPOLOGY = 1, MEMBRANE = 2
};

void exportApi(py::module &api) {
    using namespace pybind11::literals;

    exportLoopApi(api);
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
    simulation.def(py::init<std::string, ctx>())
            .def_property_readonly("single_precision", &sim::singlePrecision)
            .def_property_readonly("double_precision", &sim::doublePrecision)
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
            .def("get_selected_kernel_type", &getSelectedKernelType)
            .def("get_particle_positions", &sim::getParticlePositions)
            .def("kernel_supports_topologies", &sim::kernelSupportsTopologies)
            .def("create_topology_particle", [](sim &self, const std::string& type, readdy::Vec3 pos) {
                auto particle = self.createTopologyParticle(type, pos);
                return rpy::ReadableParticle(particle, self.context());
            }, "type"_a, "position"_a)
            .def("get_particles_for_topology", &sim::getParticlesForTopology, "topology"_a)
            .def("add_topology", [](sim &self, const std::string &name,
                                    const std::vector<rpy::ReadableParticle> &particles) {
                std::vector<readdy::model::Particle> rParticles;
                rParticles.reserve(particles.size());
                std::transform(particles.begin(), particles.end(), std::back_inserter(rParticles), [](const auto& p) {
                    return p.toParticle();
                });
                auto* top = self.addTopology(name, rParticles);
                return PyTopology(top);
            }, rvp::reference_internal, "type"_a, "particles"_a)
            .def("add_topology", [](sim &self, const std::string &name, const std::vector<std::string> &types,
                                    const py::array_t<readdy::scalar> &positions) {
                auto nParticles = positions.shape(0);
                auto nTypes = types.size();
                if(nParticles != nTypes && nTypes != 1) {
                    throw std::invalid_argument(fmt::format("the number of particles ({}) must be equal to the "
                                                                    "number of types ({})!", nParticles, nTypes));
                }
                std::vector<readdy::model::Particle> particles;
                for(std::size_t i = 0; i < nParticles; ++i) {
                    auto type =  nTypes != 1 ? types[i] : types[0];
                    particles.push_back(self.createTopologyParticle(type, readdy::Vec3(positions.at(i, 0),
                                                                                       positions.at(i, 1),
                                                                                       positions.at(i, 2))));
                }
                return PyTopology(self.addTopology(name, particles));
            }, rvp::reference_internal)
            .def_property_readonly("current_topologies", [](sim &self) {
                auto curr = self.currentTopologies();
                std::vector<PyTopology> topologies;
                topologies.reserve(curr.size());
                std::transform(curr.begin(), curr.end(), std::back_inserter(topologies), [](auto* ptr) { return PyTopology(ptr); });
                return topologies;
            })
            .def_property_readonly("current_particles", [](const sim &self) {
                std::vector<rpy::ReadableParticle> particles;
                auto currentParticles = self.currentParticles();
                particles.reserve(currentParticles.size());
                for(const auto &p : currentParticles) {
                    particles.emplace_back(p, self.context());
                }
                return particles;
            })
            .def_property_readonly("context", [](sim &self) -> const readdy::model::Context& {
                return self.context();
            })
            .def("create_loop", &sim::createLoop, py::keep_alive<0, 1>(), py::return_value_policy::reference_internal)
            .def("run", [](sim &self, const readdy::TimeStep steps, const readdy::scalar timeStep) {
                py::gil_scoped_release release;
                self.run(steps, timeStep);
            }, "n_steps"_a, "time_step"_a);
    exportObservables(api, simulation);

    // actions and evaluate observables, i.e. things needed to build a custom simulation loop [experimental]
    {
        using Action = readdy::model::actions::Action;
        using EvalObs = readdy::model::actions::EvaluateObservables;
        using MkCkpt = readdy::model::actions::MakeCheckpoint;

        auto actionsModule = api.def_submodule("actions");
        py::class_<Action>(actionsModule, "Action")
                .def("__call__", &Action::perform);

        py::class_<readdy::model::actions::top::BreakConfig>(actionsModule, "BreakConfig")
                .def(py::init<>())
                .def("add_breakable_pair", &readdy::model::actions::top::BreakConfig::addBreakablePair);

        simulation
        .def("create_action_initialize_kernel", [](sim &self) -> std::unique_ptr<Action> { return self.actions().initializeKernel(); })
        .def("create_action_euler_bd", [](sim &self, readdy::scalar timeStep) -> std::unique_ptr<Action> { return self.actions().eulerBDIntegrator(timeStep); })
        .def("create_action_calculate_forces", [](sim &self) -> std::unique_ptr<Action> { return self.actions().calculateForces();})
        .def("create_action_create_neighbor_list", [](sim &self, readdy::scalar interactionDistance) -> std::unique_ptr<Action> { return self.actions().createNeighborList(interactionDistance);})
        .def("create_action_update_neighbor_list", [](sim &self) -> std::unique_ptr<Action> { return self.actions().updateNeighborList();})
        .def("create_action_clear_neighbor_list", [](sim &self) -> std::unique_ptr<Action> { return self.actions().clearNeighborList();})
        .def("create_action_uncontrolled_approximation", [](sim &self, readdy::scalar timeStep) -> std::unique_ptr<Action> { return self.actions().uncontrolledApproximation(timeStep);})
        .def("create_action_gillespie", [](sim &self, readdy::scalar timeStep) -> std::unique_ptr<Action> { return self.actions().gillespie(timeStep);})
        .def("create_action_detailed_balance", [](sim &self, readdy::scalar timeStep) -> std::unique_ptr<Action> { return self.actions().detailedBalance(timeStep);})
        .def("create_action_evaluate_topology_reactions", [](sim &self, readdy::scalar timeStep) -> std::unique_ptr<Action> { return self.actions().evaluateTopologyReactions(timeStep);})
        .def("create_action_break_bonds", [](sim &self, readdy::scalar timeStep, const readdy::model::actions::top::BreakConfig &breakConfig) -> std::unique_ptr<Action> { return self.actions().breakBonds(timeStep, breakConfig);});

        // strictly not an action
        py::class_<EvalObs>(actionsModule, "EvaluateObservables").def("__call__", &EvalObs::perform);
        simulation.def("create_action_evaluate_observables", [](sim &self) -> std::unique_ptr<EvalObs> { return self.actions().evaluateObservables();});

        // strictly not an action
        py::class_<MkCkpt>(actionsModule, "MakeCheckpoint").def("__call__", &MkCkpt::perform);
        simulation.def("create_action_make_checkpoint", [](sim &self, const std::string &basePath, std::size_t maxNSaves) -> std::unique_ptr<MkCkpt> { return self.actions().makeCheckpoint(basePath, maxNSaves); });
    }

    struct nodelete {
        void operator()(kp* ptr) const {}
    };

    py::class_<kp, std::unique_ptr<kp, nodelete>>(api, "KernelProvider")
            .def_static("get", &kp::getInstance, rvp::reference)
            .def("load_from_dir", &kp::loadKernelsFromDirectory, "directory"_a)
            .def("available_kernels", &kp::availableKernels);

    py::class_<kern>(api, "Kernel").def("get_name", &kern::name, rvp::reference);
}
