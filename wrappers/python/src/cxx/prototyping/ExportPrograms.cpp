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
 * @file Programs.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */


#include <pybind11/pybind11.h>

#include <readdy/kernel/singlecpu/programs/SCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SCPUAddParticle.h>
#include <readdy/kernel/singlecpu/programs/SCPUEulerBDIntegrator.h>
#include <readdy/kernel/singlecpu/programs/SCPUCalculateForces.h>
#include <readdy/kernel/singlecpu/programs/SCPUUpdateNeighborList.h>
#include <readdy/kernel/singlecpu/programs/SCPUReactionImpls.h>
#include "ProgramWrap.h"

namespace py = pybind11;
namespace rpy = readdy::rpy;

using rdy_particle_t = readdy::model::Particle;

using prog_factory_t = readdy::model::programs::ProgramFactory;

using program_t = readdy::model::programs::Program;
using program_wrap_t = readdy::rpy::PyProgram;

using add_particle_t = readdy::kernel::scpu::programs::SCPUAddParticle;
using euler_integrator_t = readdy::kernel::scpu::programs::SCPUEulerBDIntegrator;
using forces_t = readdy::kernel::scpu::programs::SCPUCalculateForces;
using neighbor_list_t = readdy::kernel::scpu::programs::SCPUUpdateNeighborList;

using reactions_u_a_t = readdy::kernel::scpu::programs::reactions::SCPUUncontrolledApproximation;

void exportPrograms(py::module &proto) {
    auto f_add_particle = &prog_factory_t::createProgram<add_particle_t>;
    auto f_euler_integrator = &prog_factory_t::createProgram<euler_integrator_t>;
    auto f_forces = &prog_factory_t::createProgram<forces_t>;
    auto f_neighbor_list = &prog_factory_t::createProgram<neighbor_list_t>;
    auto f_reactions_uncontrolled_approximation = &prog_factory_t::createProgram<reactions_u_a_t>;

    using scpu_kernel_t = readdy::kernel::scpu::SCPUKernel;

    py::class_<prog_factory_t>(proto, "ProgramFactory")
            .def("create_add_particles", f_add_particle)
            .def("create_euler_integrator", f_euler_integrator)
            .def("create_update_forces", f_forces)
            .def("create_update_neighbor_list", f_neighbor_list)
            .def("create_reactions_uncontrolled_approximation", f_reactions_uncontrolled_approximation);

    py::class_ <program_t, program_wrap_t> program(proto, "Program");
    program
            .def(py::init<>())
            .def("execute", &program_t::execute);

    py::class_<add_particle_t>(proto, "AddParticle", program)
            .def(py::init<scpu_kernel_t *>())
            .def("execute", &add_particle_t::execute)
            .def("add_particle", &add_particle_t::addParticle)
            .def("set_particles", &add_particle_t::setParticles);

    py::class_<euler_integrator_t>(proto, "EulerBDIntegrator", program)
            .def(py::init<scpu_kernel_t *>())
            .def("execute", &euler_integrator_t::execute);

    py::class_<forces_t>(proto, "CalculateForces", program)
            .def(py::init<scpu_kernel_t *>())
            .def("execute", &forces_t::execute);

    py::class_<neighbor_list_t>(proto, "UpdateNeighborList", program)
            .def(py::init<scpu_kernel_t *>())
            .def("execute", &neighbor_list_t::execute);

    /**
     *
                 [](reactions_u_a_t &self, const std::string name, bpy::object handler) {
                     self.registerReactionScheme_11(name, py_fun_11_t(handler));
                 })
     */

    py::class_<reactions_u_a_t>(proto, "ReactionsUncontrolledApproximation", program)
            .def(py::init<scpu_kernel_t *>())
            .def("execute", &reactions_u_a_t::execute)
            .def("register_reaction_scheme_11", &reactions_u_a_t::registerReactionScheme_11)
            .def("register_reaction_scheme_12", &reactions_u_a_t::registerReactionScheme_12)
            .def("register_reaction_scheme_21", &reactions_u_a_t::registerReactionScheme_21)
            .def("register_reaction_scheme_22", &reactions_u_a_t::registerReactionScheme_22);
}