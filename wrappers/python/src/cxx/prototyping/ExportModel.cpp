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

using rdy_ctx_t = readdy::model::KernelContext;
using rdy_particle_t = readdy::model::Particle;

using rdy_scpu_model_t = readdy::kernel::scpu::SCPUStateModel;
using rdy_scpu_model_wrap_t = readdy::rpy::Model;

using rdy_scpu_nl_t = readdy::kernel::scpu::model::SCPUNeighborList;
using rdy_scpu_nl_box_t = readdy::kernel::scpu::model::Box;
using rdy_scpu_pd_t = readdy::kernel::scpu::model::SCPUParticleData;
using rdy_pd_entry = readdy::kernel::scpu::model::Entry;

using rdy_pot_1 = readdy::model::potentials::PotentialOrder1;
using rdy_pot_2 = readdy::model::potentials::PotentialOrder2;

void exportModelClasses(py::module &proto) {

    using rdy_uint = unsigned int;

    py::class_<rdy_particle_t>(proto, "Particle")
            .def(py::init<double, double, double, rdy_uint>())
            .def_property_readonly("pos", [](rdy_particle_t &self) { return self.getPos(); })
            .def_property_readonly("type", &rdy_particle_t::getType)
            .def_property_readonly("id", &rdy_particle_t::getId, rvp::reference_internal)
            .def(py::self == py::self)
            .def(py::self != py::self);

    py::class_<rdy_ctx_t>(proto, "Context")
            .def_property("kbt", &rdy_ctx_t::getKBT, &rdy_ctx_t::setKBT)
            .def("get_box_size", [](rdy_ctx_t &self) { return readdy::model::Vec3(self.getBoxSize()); })
            .def("set_box_size", [](rdy_ctx_t &self, readdy::model::Vec3 vec) { self.setBoxSize(vec[0], vec[1], vec[2]); })
            .def_property("periodic_boundary", &rdy_ctx_t::getPeriodicBoundary,
                          [](rdy_ctx_t &self, std::array<bool, 3> periodic) {
                              self.setPeriodicBoundary(periodic[0], periodic[1], periodic[2]);
                          })
            .def("register_particle_type", &rdy_ctx_t::registerParticleType)
            .def("get_diffusion_constant",
                 (double (rdy_ctx_t::*)(const std::string &) const) &rdy_ctx_t::getDiffusionConstant)
            .def("get_fix_position_fun", &rdy_ctx_t::getFixPositionFun)
            .def("get_shortest_difference_fun", &rdy_ctx_t::getShortestDifferenceFun)
            .def("get_dist_squared_fun", &rdy_ctx_t::getDistSquaredFun)
            .def("get_particle_radius",
                 (double (rdy_ctx_t::*)(const std::string &) const) &rdy_ctx_t::getParticleRadius)
            .def("register_conversion_reaction",
                 [](rdy_ctx_t &self, readdy::model::reactions::Conversion* r) -> const short {
                     return self.registerExternalReaction(r);
                 }, rvp::reference_internal)
            .def("register_enzymatic_reaction",
                 [](rdy_ctx_t &self, readdy::model::reactions::Enzymatic* r) -> const short {
                     return self.registerExternalReaction(r);
                 }, rvp::reference_internal)
            .def("register_fission_reaction",
                 [](rdy_ctx_t &self, readdy::model::reactions::Fission* r) -> const short {
                     return self.registerExternalReaction(r);
                 }, rvp::reference_internal)
            .def("register_fusion_reaction",
                 [](rdy_ctx_t &self, readdy::model::reactions::Fusion *r) -> const short {
                     return self.registerExternalReaction(r);
                 }, rvp::reference_internal)
            .def("register_decay_reaction",
                 [](rdy_ctx_t &self, readdy::model::reactions::Decay *r) -> const short {
                     return self.registerExternalReaction(r);
                 }, rvp::reference_internal)
            .def("register_potential_order_1",
                 [](rdy_ctx_t &self, rdy_pot_1 &pot) -> const short {
                     return self.registerExternalPotential(&pot);
                 }, rvp::reference_internal)
            .def("register_potential_order_2",
                 [](rdy_ctx_t &self, rdy_pot_2 *p) -> const short {
                     return self.registerExternalPotential(p);
                 }, rvp::reference_internal)
            .def("get_particle_type_id", &rdy_ctx_t::getParticleTypeID)
            .def("configure", &rdy_ctx_t::configure, py::arg("debugOutput") = false);

    py::class_ <rdy_scpu_model_t, rdy_scpu_model_wrap_t> model(proto, "Model");
    model
            .def(py::init<rdy_ctx_t *, readdy::model::top::TopologyActionFactory *>())
            .def("remove_particle", &rdy_scpu_model_t::removeParticle)
            .def("get_particle_positions", &rdy_scpu_model_t::getParticlePositions)
            .def("get_energy", &rdy_scpu_model_t::getEnergy)
            .def("increase_energy", &rdy_scpu_model_t::increaseEnergy)
            .def("get_particle_data", [](rdy_scpu_model_t& self) -> const rdy_scpu_pd_t& {
                return *self.getParticleData();
            }, rvp::reference_internal)
            .def("get_neighbor_list", &rdy_scpu_model_t::getNeighborList, rvp::reference_internal)
            .def("get_particles", &rdy_scpu_model_t::getParticles);

    py::class_<rdy_scpu_nl_t>(proto, "NeighborList")
            .def(py::init<rdy_ctx_t *>())
            .def("create", &rdy_scpu_nl_t::create)
            .def("setup_neighboring_boxes", &rdy_scpu_nl_t::setupNeighboringBoxes)
            .def("setup_boxes", &rdy_scpu_nl_t::setupBoxes)
            .def("fill_boxes", &rdy_scpu_nl_t::fillBoxes);
    py::class_<rdy_scpu_nl_box_t>(proto, "NeighborListBox")
            .def(py::init<long, long, long, long>())
            .def("add_neighbor", &rdy_scpu_nl_box_t::addNeighbor)
            .def_readonly("i", &rdy_scpu_nl_box_t::i)
            .def_readonly("j", &rdy_scpu_nl_box_t::j)
            .def_readonly("k", &rdy_scpu_nl_box_t::k)
            .def_readonly("id", &rdy_scpu_nl_box_t::id)
            .def_readwrite("particle_indices", &rdy_scpu_nl_box_t::particleIndices)
            .def_readwrite("neighboring_boxes", &rdy_scpu_nl_box_t::neighbors);

    py::class_<rdy_pd_entry>(proto, "ParticleDataEntry")
            .def("is_deactivated", &rdy_pd_entry::is_deactivated)
            .def("position", &rdy_pd_entry::position)
            .def_readwrite("force", &rdy_pd_entry::force)
            .def_readwrite("pos", &rdy_pd_entry::pos)
            .def_readwrite("id", &rdy_pd_entry::id)
            .def_readwrite("type", &rdy_pd_entry::type);

    py::class_<rdy_scpu_pd_t>(proto, "ParticleData")
            .def("size", &rdy_scpu_pd_t::size)
            .def("clear", &rdy_scpu_pd_t::clear)
            .def("add_particle", &rdy_scpu_pd_t::addParticle)
            .def("add_particles", &rdy_scpu_pd_t::addParticles)
            .def("remove_particle", (void (rdy_scpu_pd_t::*)(const rdy_particle_t &)) &rdy_scpu_pd_t::removeParticle)
            .def("remove_particle", (void (rdy_scpu_pd_t::*)(const std::size_t)) &rdy_scpu_pd_t::removeParticle)
            .def_property_readonly("entries", [](rdy_scpu_pd_t &self) {
                return py::make_iterator(self.begin(), self.end());
            })
            .def("__getitem__", [](rdy_scpu_pd_t &self, const unsigned int i) { return self.getParticle(i); });
}