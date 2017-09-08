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
 * @file Model.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 08.08.16
 */

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <readdy/model/Particle.h>
#include <readdy/kernel/singlecpu/model/SCPUNeighborList.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>
#include "ModelWrap.h"

namespace py = pybind11;
namespace rpy = readdy::rpy;

using rvp = py::return_value_policy;

using context = readdy::model::KernelContext;
using type_registry = readdy::model::ParticleTypeRegistry;
using potential_registry = readdy::model::potentials::PotentialRegistry;
using reaction_registry = readdy::model::reactions::ReactionRegistry;
using particle_t = readdy::model::Particle;

using model_t = readdy::kernel::scpu::SCPUStateModel;
using model_wrap = readdy::rpy::Model;

using scpu_neighbor_list = readdy::kernel::scpu::model::SCPUNeighborList;
using scpu_nl_box = readdy::kernel::scpu::model::Box;
using scpu_particle_data = readdy::kernel::scpu::model::SCPUParticleData;
using scpu_pd_entry = readdy::kernel::scpu::model::Entry;

using potential_o1 = readdy::model::potentials::PotentialOrder1;
using potential_o2 = readdy::model::potentials::PotentialOrder2;

void exportModelClasses(py::module &proto) {

    using rdy_uint = unsigned int;
    using namespace pybind11::literals;

    py::class_<particle_t>(proto, "Particle")
            .def(py::init<readdy::scalar, readdy::scalar, readdy::scalar, rdy_uint>(), "x"_a, "y"_a, "z"_a, "type"_a)
            .def_property_readonly("pos", [](particle_t &self) { return self.getPos(); })
            .def_property_readonly("type", &particle_t::getType)
            .def_property_readonly("id", &particle_t::getId, rvp::reference_internal)
            .def(py::self == py::self)
            .def(py::self != py::self);

    py::class_<type_registry>(proto, "ParticleTypeRegistry")
            .def("add", [](type_registry& self, const std::string& name, readdy::scalar D, readdy::scalar r,
                           readdy::model::particle_flavor flavor) {
                self.add(name, D, r, flavor);
            }, "name"_a, "diffusion_constant"_a, "radius"_a, "flavor"_a = 0)
            .def("diffusion_constant_of", [](const type_registry& self, const std::string& type) {
                return self.diffusion_constant_of(type);
            }, "type"_a)
            .def("radius_of", [](const type_registry& self, const std::string& type) {
                return self.radius_of(type);
            }, "type"_a)
            .def("id_of", &type_registry::id_of, "name"_a);

    py::class_<potential_registry>(proto, "PotentialRegistry")
            .def("add_external_order1", [](potential_registry& self, potential_o1& pot) {
                return self.add_external(&pot);
            })
            .def("add_external_order2", [](potential_registry& self, potential_o2& pot) {
                return self.add_external(&pot);
            });

    py::class_<context>(proto, "Context")
            .def_property("kbt", [](const context &self) {return self.kBT(); }, [](context &self, readdy::scalar kbt) {
                self.kBT() = kbt;
            })
            .def("get_box_size", [](context &self) { return readdy::model::Vec3(self.boxSize()); })
            .def("set_box_size", [](context &self, readdy::model::Vec3 vec) { self.boxSize() = vec.data; })
            .def_property("periodic_boundary", [](const context &self) {return self.periodicBoundaryConditions();},
                          [](context &self, context::PeriodicBoundaryConditions periodic) {
                              self.periodicBoundaryConditions() = periodic;
                          })
            .def("particle_types", [](context& self) -> type_registry& {
                return self.particle_types();
            }, rvp::reference_internal)
            .def("potentials", [](context& self) -> potential_registry& {
                return self.potentials();
            }, rvp::reference_internal)
            .def("get_fix_position_fun", &context::fixPositionFun)
            .def("get_shortest_difference_fun", &context::shortestDifferenceFun)
            .def("get_dist_squared_fun", &context::distSquaredFun)
            .def("register_conversion_reaction",
                 [](context &self, readdy::model::reactions::Conversion* r) -> const short {
                     return self.reactions().add_external(r);
                 }, rvp::reference_internal)
            .def("register_enzymatic_reaction",
                 [](context &self, readdy::model::reactions::Enzymatic* r) -> const short {
                     return self.reactions().add_external(r);
                 }, rvp::reference_internal)
            .def("register_fission_reaction",
                 [](context &self, readdy::model::reactions::Fission* r) -> const short {
                     return self.reactions().add_external(r);
                 }, rvp::reference_internal)
            .def("register_fusion_reaction",
                 [](context &self, readdy::model::reactions::Fusion *r) -> const short {
                     return self.reactions().add_external(r);
                 }, rvp::reference_internal)
            .def("register_decay_reaction",
                 [](context &self, readdy::model::reactions::Decay *r) -> const short {
                     return self.reactions().add_external(r);
                 }, rvp::reference_internal)
            .def("configure", &context::configure, py::arg("debugOutput") = false);

    py::class_ <model_t, model_wrap> model(proto, "Model");
    model
            .def(py::init<context *, readdy::model::top::TopologyActionFactory *>())
            .def("remove_particle", &model_t::removeParticle)
            .def("get_particle_positions", &model_t::getParticlePositions)
            .def("get_energy", &model_t::getEnergy)
            .def("increase_energy", &model_t::increaseEnergy)
            .def("get_particle_data", [](model_t& self) -> const scpu_particle_data& {
                return *self.getParticleData();
            }, rvp::reference_internal)
            .def("get_neighbor_list", &model_t::getNeighborList, rvp::reference_internal)
            .def("get_particles", &model_t::getParticles);

    py::class_<scpu_neighbor_list>(proto, "NeighborList")
            .def(py::init<context *>())
            .def("create", &scpu_neighbor_list::create)
            .def("setup_neighboring_boxes", &scpu_neighbor_list::setupNeighboringBoxes)
            .def("setup_boxes", &scpu_neighbor_list::setupBoxes)
            .def("fill_boxes", &scpu_neighbor_list::fillBoxes);
    py::class_<scpu_nl_box>(proto, "NeighborListBox")
            .def(py::init<long, long, long, long>())
            .def("add_neighbor", &scpu_nl_box::addNeighbor)
            .def_readonly("i", &scpu_nl_box::i)
            .def_readonly("j", &scpu_nl_box::j)
            .def_readonly("k", &scpu_nl_box::k)
            .def_readonly("id", &scpu_nl_box::id)
            .def_readwrite("particle_indices", &scpu_nl_box::particleIndices)
            .def_readwrite("neighboring_boxes", &scpu_nl_box::neighbors);

    py::class_<scpu_pd_entry>(proto, "ParticleDataEntry")
            .def("is_deactivated", &scpu_pd_entry::is_deactivated)
            .def("position", &scpu_pd_entry::position)
            .def_readwrite("force", &scpu_pd_entry::force)
            .def_readwrite("pos", &scpu_pd_entry::pos)
            .def_readwrite("id", &scpu_pd_entry::id)
            .def_readwrite("type", &scpu_pd_entry::type);

    py::class_<scpu_particle_data>(proto, "ParticleData")
            .def("size", &scpu_particle_data::size)
            .def("clear", &scpu_particle_data::clear)
            .def("add_particle", &scpu_particle_data::addParticle)
            .def("add_particles", &scpu_particle_data::addParticles)
            .def("remove_particle", (void (scpu_particle_data::*)(const particle_t &)) &scpu_particle_data::removeParticle)
            .def("remove_particle", (void (scpu_particle_data::*)(const std::size_t)) &scpu_particle_data::removeParticle)
            .def_property_readonly("entries", [](scpu_particle_data &self) {
                return py::make_iterator(self.begin(), self.end());
            })
            .def("__getitem__", [](scpu_particle_data &self, const unsigned int i) { return self.getParticle(i); });
}