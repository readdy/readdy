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
 * @file ExportSchemeApi.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.08.16
 */

#ifndef READDY_MAIN_EXPORTSCHEMEAPI_H
#define READDY_MAIN_EXPORTSCHEMEAPI_H

#include <pybind11/pybind11.h>
#include <readdy/api/SimulationLoop.h>
#include <readdy/model/actions/UserDefinedAction.h>
#include "PyFunction.h"

using UserAction = readdy::model::actions::UserDefinedAction;


void exportSchemeApi(pybind11::module &module) {
    namespace py = pybind11;
    using namespace py::literals;
    using Loop = readdy::api::SimulationLoop;

    py::class_<UserAction, std::shared_ptr<UserAction>> userAction (module, "UserDefinedAction");

    py::class_<Loop>(module, "SimulationLoop")
            .def_property("progress_callback", [](const Loop& self) { return self.progressCallback(); },
                          [](Loop &self, const std::function<void(readdy::time_step_type)> &fun) {
                self.progressCallback() = fun;
            })
            .def_property("progress_output_stride", [](const Loop& self) { return self.progressOutputStride(); },
                          [](Loop &self, std::size_t stride) { self.progressOutputStride() = stride; })
            .def("run", [](Loop& self, const readdy::time_step_type steps) {
                py::gil_scoped_release release;
                self.run(steps);
            }, "n_steps"_a)
            .def("run_with_criterion", [](Loop& self, pybind11::object continuingCriterion) {
                py::gil_scoped_release release;
                auto pyFun = readdy::rpy::PyFunction<bool(const readdy::time_step_type current)>(continuingCriterion);
                self.run(pyFun);
            }, "continuing_criterion"_a, py::keep_alive<0, 1>())
            .def("run_initialize", &Loop::runInitialize)
            .def("run_initialize_neighbor_list", &Loop::runInitializeNeighborList)
            .def("run_update_neighbor_list", &Loop::runUpdateNeighborList)
            .def("run_clear_neighbor_list", &Loop::runClearNeighborList)
            .def("run_forces", &Loop::runForces)
            .def("run_evaluate_observables", &Loop::runEvaluateObservables)
            .def("run_integrator", &Loop::runIntegrator)
            .def("run_reactions", &Loop::runReactions)
            .def("run_topology_reactions", &Loop::runTopologyReactions)
            .def("use_integrator", [](Loop &self, std::string name) {
                self.useIntegrator(name, self.timeStep());
            })
            .def("use_integrator", [](Loop &self, std::shared_ptr<UserAction> integrator) {
                integrator->kernel() = self.kernel();
                self.integrator() = integrator;
            }, py::keep_alive<1, 2>())
            .def("evaluate_forces", &Loop::evaluateForces, "evaluate"_a)
            .def("use_reaction_scheduler", [](Loop &self, std::string name) {
                     return self.useReactionScheduler(name, self.timeStep());
                 }, "reaction_scheduler_name"_a)
            .def("use_reaction_scheduler", [](Loop &self, std::shared_ptr<UserAction> reactionScheduler) {
                reactionScheduler->kernel() = self.kernel();
                self.reactionScheduler() = reactionScheduler;
            }, py::keep_alive<1, 2>())
            .def("write_config_to_file", &Loop::writeConfigToFile, py::return_value_policy::reference_internal, "file"_a)
            .def("evaluate_topology_reactions", [](Loop &self, bool evaluate, py::object timeStep) {
                self.evaluateTopologyReactions(evaluate, timeStep.is_none() ? self.timeStep() : timeStep.cast<readdy::scalar>());
            }, "evaluate"_a, "timeStep"_a = py::none())
            .def("evaluate_observables", &Loop::evaluateObservables, "evaluate"_a)
            .def_property("skin_size", [](const Loop &self) { return self.skinSize(); },
                          [](Loop &self, readdy::scalar skin) { self.skinSize() = skin; })
            .def("validate", &Loop::validate);
}

#endif //READDY_MAIN_EXPORTSCHEMEAPI_H
