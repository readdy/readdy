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
#include <readdy/api/SimulationScheme.h>
#include "PyFunction.h"


template<typename SchemeType>
void exportSchemeApi(pybind11::module &module, std::string schemeName) {
    namespace py = pybind11;
    using namespace py::literals;
    using conf = readdy::api::SchemeConfigurator<SchemeType>;
    py::class_<SchemeType>(module, schemeName.c_str())
            .def("run", [](SchemeType& self, const readdy::time_step_type steps) {
                py::gil_scoped_release release;
                self.run(steps);
            }, "n_steps"_a)
            .def("run_with_criterion", [](SchemeType& self, pybind11::object continuingCriterion) {
                py::gil_scoped_release release;
                auto pyFun = readdy::rpy::PyFunction<bool(const readdy::time_step_type current)>(continuingCriterion);
                self.run(std::move(pyFun));
            }, "continuing_criterion"_a);
    std::string configuratorName =  schemeName + "Configurator";
    py::class_<conf>(module, configuratorName.c_str())
            .def("with_integrator",
                 [](conf &self, std::string name) -> conf & { return self.withIntegrator(name); },
                 py::return_value_policy::reference_internal, "integrator_name"_a)
            .def("include_forces", &conf::includeForces, py::return_value_policy::reference_internal,
                 "do_include"_a = true)
            .def("with_reaction_scheduler",
                 [](conf &self, std::string name) -> conf & {
                     return (self.withReactionScheduler(name));
                 },
                 py::return_value_policy::reference_internal, "reaction_scheduler_name"_a)
            .def("write_config_to_file", &conf::writeConfigToFile, py::return_value_policy::reference_internal, "file"_a)
            .def("evaluate_observables", &conf::evaluateObservables, py::return_value_policy::reference_internal,
                 "do_evaluate"_a = true)
            .def("with_skin_size", &conf::withSkinSize, py::return_value_policy::reference_internal, "skin_size"_a = -1)
            .def("configure", &conf::configure, "time_step"_a)
            .def("configure_and_run", [](conf& self, double dt, const readdy::time_step_type steps) {
                py::gil_scoped_release release;
                self.configureAndRun(dt, steps);
            }, "time_step"_a, "n_time_steps"_a);
}

template<>
void exportSchemeApi<readdy::api::AdvancedScheme>(pybind11::module &module, std::string schemeName) {
    namespace py = pybind11;
    using namespace py::literals;
    using scheme_t = readdy::api::AdvancedScheme;
    using conf = readdy::api::SchemeConfigurator<scheme_t>;
    py::class_<scheme_t>(module, schemeName.c_str())
            .def("run", [](scheme_t& self, const readdy::time_step_type steps) {
                py::gil_scoped_release release;
                self.run(steps);
            }, "n_steps"_a)
            .def("run_with_criterion", [](scheme_t& self, pybind11::object continuingCriterion) {
                py::gil_scoped_release release;
                auto pyFun = readdy::rpy::PyFunction<bool(const readdy::time_step_type current)>(continuingCriterion);
                self.run(std::move(pyFun));
            }, "continue_criterion"_a);
    std::string configuratorName =  schemeName + "Configurator";
    py::class_<conf>(module, configuratorName.c_str())
            .def("with_integrator",
                 [](conf &self, std::string name) -> conf & { return self.withIntegrator(name); },
                 py::return_value_policy::reference_internal, "integrator_name"_a)
            .def("include_forces", &conf::includeForces, py::return_value_policy::reference_internal, "do_include"_a = true)
            .def("include_compartments", &conf::includeCompartments, py::return_value_policy::reference_internal, "do_include"_a = true)
            .def("with_reaction_scheduler",
                 [](conf &self, std::string name) -> conf & {
                     return (self.withReactionScheduler(name));
                 },
                 py::return_value_policy::reference_internal, "reaction_scheduler_name"_a)
            .def("evaluate_observables", &conf::evaluateObservables, py::return_value_policy::reference_internal, "do_evaluate"_a = true)
            .def("with_skin_size", &conf::withSkinSize, py::return_value_policy::reference_internal, "skin_size"_a = -1)
            .def("configure", &conf::configure, "time_step"_a)
            .def("configure_and_run", [](conf& self, double dt, const readdy::time_step_type steps) {
                py::gil_scoped_release release;
                self.configureAndRun(dt, steps);
            }, "time_step"_a, "n_time_steps"_a);
};

#endif //READDY_MAIN_EXPORTSCHEMEAPI_H
