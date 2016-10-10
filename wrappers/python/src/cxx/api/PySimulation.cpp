#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <boost/uuid/uuid_io.hpp>
#include <readdy/Simulation.h>
#include <readdy/plugin/KernelProvider.h>
#include "ExportSchemeApi.h"
#include "PyPotential.h"
#include "PyFunction.h"

namespace bpy = pybind11;

using rvp = bpy::return_value_policy;
using sim = readdy::Simulation;
using kp = readdy::plugin::KernelProvider;
using vec = readdy::model::Vec3;
using pot2 = readdy::py::PotentialOrder2Wrapper;
using model = readdy::model::KernelStateModel;
using ctx = readdy::model::KernelContext;
using kern = readdy::model::Kernel;
using uuid = boost::uuids::uuid;

struct nodelete {
    template<typename T>
    void operator()(T *) {}
};

// thin wrappers
void setBoxSize(sim &self, const vec &size) { /* explicitly choose void(vec) signature */ self.setBoxSize(size); }

std::string getSelectedKernelType(sim &self) { /* discard const reference */ return self.getSelectedKernelType(); }

void addParticle(sim &self, const std::string &type, const vec &pos) { self.addParticle(pos[0], pos[1], pos[2], type); }

void registerPotentialOrder2(sim &self, pot2 *potential, std::string type1, std::string type2) {
    self.registerPotentialOrder2(potential, type1, type2);
}

boost::uuids::uuid
registerObservable_ParticlePositions(sim &self, unsigned int stride, pybind11::object callbackFun) {
    auto pyFun = readdy::py::PyFunction<void(readdy::model::ParticlePositionObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::ParticlePositionObservable>(std::move(pyFun), stride);
}

boost::uuids::uuid
registerObservable_RadialDistribution(sim &self, unsigned int stride, pybind11::object callbackFun,
                                      bpy::array_t<double> &binBorders, std::string typeCountFrom,
                                      std::string typeCountTo, double particleDensity) {
    auto pyFun = readdy::py::PyFunction<void(readdy::model::RadialDistributionObservable::result_t)>(callbackFun);
    const auto info = binBorders.request();
    std::vector<double> binBordersVec{};
    binBordersVec.reserve(info.shape[0]);
    const auto data = static_cast<double *>(info.ptr);
    for (auto i = 0; i < info.shape[0]; ++i) binBordersVec.push_back(data[i]);
    return self.registerObservable<readdy::model::RadialDistributionObservable>(std::move(pyFun), stride, binBordersVec,
                                                                                typeCountFrom, typeCountTo,
                                                                                particleDensity);
}

boost::uuids::uuid
registerObservable_CenterOfMass(sim &self, unsigned int stride, const pybind11::object &callbackFun,
                                std::vector<std::string> types) {
    auto pyFun = readdy::py::PyFunction<void(readdy::model::CenterOfMassObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::CenterOfMassObservable>(
            std::move(pyFun), stride, types
    );
}

boost::uuids::uuid
registerObservable_HistogramAlongAxisObservable(sim &self, unsigned int stride, const bpy::object &callbackFun,
                                                bpy::array_t<double> &binBorders, bpy::list types, unsigned int axis) {
    const auto info = binBorders.request();
    const auto sizeBorders = info.shape[0];
    auto binBordersData = static_cast<double *>(info.ptr);
    const auto sizeTypes = bpy::len(types);
    std::vector<std::string> typesVec{};
    typesVec.reserve((unsigned long) sizeTypes);
    std::vector<double> binBordersVec{};
    binBordersVec.reserve(sizeBorders);
    for (auto i = 0; i < sizeBorders; ++i) {
        binBordersVec.push_back(binBordersData[i]);
    }
    for (auto i = 0; i < sizeTypes; ++i) {
        typesVec.push_back(types[i].cast<std::string>());
    }
    auto pyFun = readdy::py::PyFunction<void(readdy::model::HistogramAlongAxisObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::HistogramAlongAxisObservable>(std::move(pyFun), stride, binBordersVec,
                                                                                typesVec, axis);

}

boost::uuids::uuid
registerObservable_NParticlesTypes(sim &self, unsigned int stride, bpy::list types, const bpy::object &callbackFun) {
    const auto sizeTypes = bpy::len(types);
    std::vector<std::string> typesVec{};
    typesVec.reserve((unsigned long) sizeTypes);
    for (auto i = 0; i < sizeTypes; ++i) {
        typesVec.push_back(types[i].cast<std::string>());
    }
    auto pyFun = readdy::py::PyFunction<void(readdy::model::NParticlesObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::NParticlesObservable>(std::move(pyFun), stride, typesVec);
}

boost::uuids::uuid registerObservable_NParticles(sim &self, unsigned int stride, const bpy::object &callbackFun) {
    auto pyFun = readdy::py::PyFunction<void(readdy::model::NParticlesObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::NParticlesObservable>(std::move(pyFun), stride);
}

// todo @chrisfroe register forces observable

// module
PYBIND11_PLUGIN (api) {

    bpy::module api("api", "ReaDDy c++-api python module");

    exportSchemeApi<readdy::api::ReaDDyScheme>(api, "ReaDDyScheme");

    bpy::class_<sim>(api, "Simulation")
            .def(bpy::init<>())
            .def_property("kbt", &sim::getKBT, &sim::setKBT)
            .def_property("periodic_boundary", &sim::getPeriodicBoundary, &sim::setPeriodicBoundary)
            .def_property("box_size", &sim::getBoxSize, &setBoxSize)
            .def("register_particle_type", &sim::registerParticleType)
            .def("add_particle", [](sim &self, const std::string &type, const vec &pos) {
                self.addParticle(pos[0], pos[1], pos[2], type);
            })
            .def("is_kernel_selected", &sim::isKernelSelected)
            .def("get_selected_kernel_type", &getSelectedKernelType)
            .def("register_potential_order_2", &registerPotentialOrder2)
            .def("register_potential_harmonic_repulsion", &sim::registerHarmonicRepulsionPotential)
            .def("register_potential_piecewise_weak_interaction",
                 &sim::registerWeakInteractionPiecewiseHarmonicPotential)
            .def("register_potential_box", &sim::registerBoxPotential)
            .def("get_particle_positions", &sim::getParticlePositions)
            .def("register_observable_particle_positions", &registerObservable_ParticlePositions)
            .def("register_observable_radial_distribution", &registerObservable_RadialDistribution)
            .def("register_observable_histogram_along_axis", &registerObservable_HistogramAlongAxisObservable)
            .def("register_observable_center_of_mass", &registerObservable_CenterOfMass)
            .def("register_observable_n_particles", &registerObservable_NParticles)
            .def("register_observable_n_particles_types", &registerObservable_NParticlesTypes)
            .def("register_reaction_conversion", &sim::registerConversionReaction, rvp::reference_internal)
            .def("register_reaction_enzymatic", &sim::registerEnzymaticReaction, rvp::reference_internal)
            .def("register_reaction_fission", &sim::registerFissionReaction, rvp::reference_internal)
            .def("register_reaction_fusion", &sim::registerFusionReaction, rvp::reference_internal)
            .def("register_reaction_decay", &sim::registerDecayReaction, rvp::reference_internal)
            .def("get_recommended_time_step", &sim::getRecommendedTimeStep)
            .def("set_kernel", &sim::setKernel)
            .def("set_time_step", &sim::setTimeStep)
            .def("run_scheme_readdy", [](sim &self, bool defaults) {
                     return std::make_unique<readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>>(
                             self.runScheme<readdy::api::ReaDDyScheme>(defaults)
                     );
                 }
            )
            .def("run", &sim::run);

    bpy::class_<kp, std::unique_ptr<kp, nodelete>>(api, "KernelProvider")
            .def_static("get", &kp::getInstance, rvp::reference)
            .def("load_from_dir", &kp::loadKernelsFromDirectory);

    bpy::class_<pot2>(api, "Pot2")
            .def(bpy::init<std::string, bpy::object, bpy::object>())
            .def("calc_energy", &pot2::calculateEnergy)
            .def("calc_force", &pot2::calculateForce);

    bpy::class_<kern>(api, "Kernel").def("get_name", &kern::getName, rvp::reference);

    return api.ptr();

}
