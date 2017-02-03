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
void exportReaDDySchemeApi(pybind11::module &module, std::string schemeName) {
    namespace py = pybind11;
    using conf = readdy::api::SchemeConfigurator<SchemeType>;
    py::class_<SchemeType>(module, schemeName.c_str())
            .def("run", [](SchemeType& self, const readdy::model::observables::time_step_type steps) {
                py::gil_scoped_release release;
                self.run(steps);
            });
    std::string configuratorName =  schemeName + "Configurator";
    py::class_<conf>(module, configuratorName.c_str())
            .def("with_integrator",
                 [](conf &self, std::string name) -> conf & { return self.withIntegrator(name); },
                 py::return_value_policy::reference_internal)
            .def("include_forces", &conf::includeForces, py::return_value_policy::reference_internal)
            .def("with_reaction_scheduler",
                 [](conf &self, std::string name) -> conf & {
                     return (self.withReactionScheduler(name));
                 },
                 py::return_value_policy::reference_internal)
            .def("evaluate_observables", &conf::evaluateObservables, py::return_value_policy::reference_internal)
            .def("configure", &conf::configure)
            .def("configure_and_run", [](conf& self, double dt, const readdy::model::observables::time_step_type steps) {
                py::gil_scoped_release release;
                self.configureAndRun(dt, steps);
            });
}

template<typename SchemeType>
void exportAdvancedSchemeApi(pybind11::module &module, std::string schemeName) {
    namespace py = pybind11;
    using conf = readdy::api::AdvancedSchemeConfigurator<SchemeType>;
    py::class_<SchemeType>(module, schemeName.c_str())
            .def("run", [](SchemeType& self, const readdy::model::observables::time_step_type steps) {
                py::gil_scoped_release release;
                self.run(steps);
            });
    std::string configuratorName =  schemeName + "Configurator";
    py::class_<conf>(module, configuratorName.c_str())
            .def("with_integrator",
                 [](conf &self, std::string name) -> conf & { return self.withIntegrator(name); },
                 py::return_value_policy::reference_internal)
            .def("include_forces", &conf::includeForces, py::return_value_policy::reference_internal)
            .def("include_compartments", &conf::includeCompartments, py::return_value_policy::reference_internal)
            .def("with_reaction_scheduler",
                 [](conf &self, std::string name) -> conf & {
                     return (self.withReactionScheduler(name));
                 },
                 py::return_value_policy::reference_internal)
            .def("evaluate_observables", &conf::evaluateObservables, py::return_value_policy::reference_internal)
            .def("configure", &conf::configure)
            .def("configure_and_run", [](conf& self, double dt, const readdy::model::observables::time_step_type steps) {
                py::gil_scoped_release release;
                self.configureAndRun(dt, steps);
            });
}


#endif //READDY_MAIN_EXPORTSCHEMEAPI_H
