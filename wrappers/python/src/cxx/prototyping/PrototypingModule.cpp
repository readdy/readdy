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

using core_kernel_t = readdy::model::Kernel;
using core_kernel_wrap_t = readdy::rpy::KernelWrap;
using core_program_factory = readdy::model::actions::ActionFactory;
using core_program_t = readdy::model::actions::Action;

void exportPrototyping(py::module& proto) {
    exportPrograms(proto);
    exportModelClasses(proto);
    exportPotentials(proto);

    py::class_<scpu_kernel_t>(proto, "SingleCPUKernel")
            .def(py::init<>())
            .def("get_kernel_state_model", [](const scpu_kernel_t &self) -> const scpu::SCPUStateModel& {return self.getSCPUKernelStateModel(); }, rvp::reference_internal)
            .def("get_kernel_context", [](const scpu_kernel_t &self) -> const readdy::model::KernelContext& { return self.getKernelContext(); }, rvp::reference_internal)
            .def("get_available_potentials", &scpu_kernel_t::getAvailablePotentials)
            .def("get_potential_factory", [](const scpu_kernel_t &self) -> const readdy::model::potentials::PotentialFactory& {return self.getPotentialFactory(); }, rvp::reference_internal)
            .def("get_reaction_factory", [](const scpu_kernel_t &self) -> const readdy::model::reactions::ReactionFactory& {return self.getReactionFactory();}, rvp::reference_internal)
            .def("get_observable_factory", [](const scpu_kernel_t &self) -> const readdy::model::observables::ObservableFactory& {return self.getObservableFactory();}, rvp::reference_internal)
            .def("get_topology_action_factory", [](const scpu_kernel_t &self) -> const readdy::model::top::TopologyActionFactory*  {return self.getTopologyActionFactory();}, rvp::reference_internal)
            .def("get_action_factory", [](const scpu_kernel_t &self) -> const readdy::model::actions::ActionFactory& {return self.getActionFactory();}, rvp::reference_internal);

}
