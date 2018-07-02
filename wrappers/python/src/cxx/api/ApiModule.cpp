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
#include <readdy/model/actions/UserDefinedAction.h>
#include "ExportObservables.h"

namespace py = pybind11;

using rvp = py::return_value_policy;
using sim = readdy::Simulation;
using kp = readdy::plugin::KernelProvider;
using vec = readdy::Vec3;
using model = readdy::model::StateModel;
using ctx = readdy::model::Context;
using kern = readdy::model::Kernel;

void exportTopologies(py::module &);
void exportSchemeApi(py::module &);
void exportKernelContext(py::module &);

std::string getSelectedKernelType(sim &self) { /* discard const reference */ return self.selectedKernelType(); }

void addParticle(sim &self, const std::string &type, const vec &pos) { self.addParticle(type, pos[0], pos[1], pos[2]); }


enum class ParticleTypeFlavor {
    NORMAL = 0, TOPOLOGY = 1, MEMBRANE = 2
};

void exportApi(py::module &api) {
    using namespace pybind11::literals;

    exportSchemeApi(api);
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
    simulation.def(py::init<std::string>())
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
            .def("create_topology_particle", &sim::createTopologyParticle, "type"_a, "position"_a)
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
            .def_property("context", [](sim &self) -> readdy::model::Context & {
                return self.context();
            }, [](sim &self, const readdy::model::Context &context) {
                self.context() = context;
            })
            .def("create_loop", &sim::createLoop, py::keep_alive<0, 1>(), py::return_value_policy::reference_internal)
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

    py::class_<kern>(api, "Kernel").def("get_name", &kern::name, rvp::reference);
}
