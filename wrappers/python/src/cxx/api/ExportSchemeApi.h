/**
 * << detailed description >>
 *
 * @file ExportSchemeApi.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.08.16
 */

#ifndef READDY_MAIN_EXPORTSCHEMEAPI_H
#define READDY_MAIN_EXPORTSCHEMEAPI_H

#include <pybind11/pybind11.h>
#include <readdy/SimulationScheme.h>


template<typename SchemeType>
void exportSchemeApi(pybind11::module &module, std::string schemeName) {
    namespace py = pybind11;
    using conf = readdy::api::SchemeConfigurator<SchemeType>;
    py::class_<SchemeType>(module, schemeName.c_str())
            .def("run", [](SchemeType& self, const readdy::model::observables::time_step_type steps) {
                py::gil_scoped_release release;
                self.run(steps);
            });
    std::string configuratorName = "SchemeConfigurator" + schemeName;
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
            .def("configure_and_run", [](conf& self, const readdy::model::observables::time_step_type steps) {
                py::gil_scoped_release release;
                self.configureAndRun(steps);
            });
}

#endif //READDY_MAIN_EXPORTSCHEMEAPI_H
