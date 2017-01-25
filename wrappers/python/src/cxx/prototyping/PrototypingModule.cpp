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
 * @file PrototypingModule.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.08.16
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include "KernelWrap.h"
#include "../api/PyFunction.h"

namespace py = pybind11;

using rvp = py::return_value_policy;

void exportPrograms(py::module &);

void exportModelClasses(py::module &);

void exportPotentials(py::module &);

namespace scpu = readdy::kernel::scpu;

using rdy_scpu_model_t = scpu::SCPUStateModel;
using scpu_kernel_t = scpu::SCPUKernel;

using scpu_kernel_wrap_t = scpu_kernel_t; // todo: do i need readdy::py::SingleCPUKernelWrap here?

using core_kernel_t = readdy::model::Kernel;
using core_kernel_wrap_t = readdy::rpy::KernelWrap;
using core_program_factory = readdy::model::actions::ActionFactory;
using core_program_t = readdy::model::actions::Action;


// module
PYBIND11_PLUGIN (prototyping) {

    py::module proto("prototyping", "ReaDDy prototyping python module");

    exportPrograms(proto);
    exportModelClasses(proto);
    exportPotentials(proto);

    py::class_<scpu_kernel_wrap_t>(proto, "SingleCPUKernel")
            .def(py::init<>())
            .def("get_kernel_state_model", &scpu_kernel_wrap_t::getKernelStateModel, rvp::reference_internal)
            .def("get_kernel_context", &scpu_kernel_wrap_t::getKernelContext, rvp::reference_internal)
            .def("get_available_potentials", &scpu_kernel_wrap_t::getAvailablePotentials)
            .def("get_potential_factory", &scpu_kernel_wrap_t::getPotentialFactory, rvp::reference_internal)
            .def("get_reaction_factory", &scpu_kernel_wrap_t::getReactionFactory, rvp::reference_internal)
            .def("get_observable_factory", &scpu_kernel_wrap_t::getObservableFactory, rvp::reference_internal)
            .def("get_program_factory", &scpu_kernel_wrap_t::getActionFactory, rvp::reference_internal);

    return proto.ptr();
}