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
#include <readdy/model/reactions/Conversion.h>
#include <readdy/model/reactions/Enzymatic.h>
#include <readdy/model/reactions/Fission.h>
#include <readdy/model/reactions/Fusion.h>
#include <readdy/model/reactions/Decay.h>

namespace py = pybind11;

using rvp = py::return_value_policy;

using context = readdy::model::Context;
using particle_t = readdy::model::Particle;

using model_t = readdy::kernel::scpu::SCPUStateModel;

using scpu_neighbor_list = readdy::kernel::scpu::model::CellLinkedList;
using scpu_particle_data = readdy::kernel::scpu::model::SCPUParticleData;
using scpu_pd_entry = readdy::kernel::scpu::model::Entry;

using potential_o1 = readdy::model::potentials::PotentialOrder1;
using potential_o2 = readdy::model::potentials::PotentialOrder2;

void exportModelClasses(py::module &proto) {

    using rdy_uint = unsigned int;
    using namespace pybind11::literals;

    py::class_<particle_t>(proto, "Particle")
            .def(py::init<readdy::scalar, readdy::scalar, readdy::scalar, readdy::particle_type_type>(), "x"_a, "y"_a, "z"_a, "type"_a)
            .def_property_readonly("pos", [](particle_t &self) { return self.getPos(); })
            .def_property_readonly("type", &particle_t::getType)
            .def_property_readonly("id", &particle_t::getId, rvp::reference_internal)
            .def(py::self == py::self)
            .def(py::self != py::self);

    py::class_ <model_t> model(proto, "Model");
    model
            .def(py::init<const context&, readdy::model::top::TopologyActionFactory *>())
            .def("remove_particle", &model_t::removeParticle)
            .def("get_particle_positions", &model_t::getParticlePositions)
            .def("get_energy", [](const model_t &self) {return self.energy();})
            .def("increase_energy", &model_t::increaseEnergy)
            .def("get_particle_data", [](model_t& self) -> const scpu_particle_data& {
                return *self.getParticleData();
            }, rvp::reference_internal)
            .def("get_particles", &model_t::getParticles);

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