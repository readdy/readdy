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
    namespace bpy = pybind11;
    using conf = readdy::api::SchemeConfigurator<SchemeType>;
    bpy::class_<SchemeType>(module, schemeName.c_str()).def("run", &SchemeType::run);
    std::string configuratorName = "SchemeConfigurator" + schemeName;
    bpy::class_<conf, std::shared_ptr<conf>>(module, configuratorName.c_str())
            .def("with_integrator",
                 [](conf &self, std::string name) -> conf & { return self.withIntegrator(name); },
                 bpy::return_value_policy::reference_internal)
            .def("include_forces", &conf::includeForces, bpy::return_value_policy::reference_internal)
            .def("with_reaction_scheduler",
                 [](conf &self, std::string name) -> conf & {
                     return (self.withReactionScheduler(name));
                 },
                 bpy::return_value_policy::reference_internal)
            .def("evaluate_observables", &conf::evaluateObservables, bpy::return_value_policy::reference_internal)
            .def("configure", &conf::configure)
            .def("configure_and_run", &conf::configureAndRun);
}

#endif //READDY_MAIN_EXPORTSCHEMEAPI_H
