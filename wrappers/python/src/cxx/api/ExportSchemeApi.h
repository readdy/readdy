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

#include <boost/python.hpp>
#include <readdy/SimulationScheme.h>


template<typename SchemeType>
void exportSchemeApi(std::string schemeName) {
    namespace bpy = boost::python;
    using configurator = readdy::api::SchemeConfigurator<SchemeType>;
    bpy::class_<SchemeType, boost::noncopyable, SchemeType*>(schemeName.c_str(), bpy::no_init).def("run", &SchemeType::run);
    std::string configuratorName = "SchemeConfigurator"+schemeName;
    bpy::class_<configurator, boost::noncopyable, std::shared_ptr<configurator>>(configuratorName.c_str(), bpy::no_init)
            .def("with_integrator", +[](configurator& self, std::string name) -> configurator& {return self.withIntegrator(name);}, bpy::return_internal_reference<>())
            .def("include_forces", &configurator::includeForces, bpy::return_internal_reference<>())
            .def("with_reaction_scheduler", +[](configurator& self, std::string name) -> configurator& {return (self.withReactionScheduler(name));}, bpy::return_internal_reference<>())
            .def("evaluate_observables", &configurator::evaluateObservables, bpy::return_internal_reference<>())
            .def("configure", readdy::py::adapt_unique(&configurator::configure))
            .def("configure_and_run", &configurator::configureAndRun);
}
#endif //READDY_MAIN_EXPORTSCHEMEAPI_H
