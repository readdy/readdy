
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include "../PyConverters.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <readdy/Simulation.h>
#include <readdy/plugin/KernelProvider.h>
#include "../PyPotential.h"
#include "../PyFunction.h"

namespace bpy = boost::python;
using sim = readdy::Simulation;
using kp = readdy::plugin::KernelProvider;
using vec = readdy::model::Vec3;
using pot2 = readdy::py::PotentialOrder2Wrapper;
using model = readdy::model::KernelStateModel;
using ctx = readdy::model::KernelContext;
using kern = readdy::model::Kernel;
using uuid = boost::uuids::uuid;

boost::python::list getPeriodicBoundarySimulationWrapper(const sim &simulation) {
    auto p = simulation.getPeriodicBoundary();
    bpy::list result;
    result.append<bool>(p[0]);
    result.append<bool>(p[1]);
    result.append<bool>(p[2]);
    return result;
}

void setPeriodicBoundarySimulationWrapper(sim &self, bpy::list list) {
    self.setPeriodicBoundary({bpy::extract<bool>(list[0]), bpy::extract<bool>(list[1]), bpy::extract<bool>(list[2])});
}

// thin wrappers
void setBoxSize(sim &self, const vec& size) { /* explicitly choose void(vec) signature */ self.setBoxSize(size); }
std::string getSelectedKernelType(sim &self) { /* discard const reference */ return self.getSelectedKernelType(); }
void addParticle(sim& self, const std::string& type, const vec& pos) { self.addParticle(pos[0], pos[1], pos[2], type); }
void registerPotentialOrder2(sim& self, pot2& potential, std::string type1, std::string type2) {
    std::unique_ptr<pot2> ptr {potential.replicate()};
    self.registerPotentialOrder2(ptr.get(), type1, type2);
}

double pyVec3Bracket(vec& self, const unsigned int i) {return self[i];}

boost::uuids::uuid registerObservable_ParticlePositions(sim& self, unsigned int stride, const boost::python::object &callbackFun) {
    auto pyFun = readdy::py::PyFunction<void(readdy::model::ParticlePositionObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::ParticlePositionObservable>(std::move(pyFun), stride);
}

boost::uuids::uuid registerObservable_RadialDistribution(sim& self, unsigned int stride, const boost::python::object &callbackFun, boost::python::numeric::array& binBorders, std::string typeCountFrom, std::string typeCountTo, double particleDensity) {
    auto pyFun = readdy::py::PyFunction<void(readdy::model::RadialDistributionObservable::result_t)>(callbackFun);
    const auto size = boost::python::len(binBorders);
    std::vector<double> binBordersVec {};
    binBordersVec.reserve((unsigned long) size);
    for(auto i = 0; i < size; ++i) binBordersVec.push_back(boost::python::extract<double>(binBorders[i]));
    return self.registerObservable<readdy::model::RadialDistributionObservable>(std::move(pyFun), stride, binBordersVec, typeCountFrom, typeCountTo, particleDensity);
}

boost::uuids::uuid registerObservable_CenterOfMass(sim& self, unsigned int stride, const boost::python::object &callbackFun, boost::python::list types) {
    auto pyFun = readdy::py::PyFunction<void(readdy::model::CenterOfMassObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::CenterOfMassObservable>(
            std::move(pyFun), stride, readdy::py::sequence_to_vector<std::string>(types)
    );
}

boost::uuids::uuid registerObservable_HistogramAlongAxisObservable(sim& self, unsigned int stride, const bpy::object& callbackFun, bpy::numeric::array& binBorders, bpy::list types, unsigned int axis) {
    const auto sizeBorders = bpy::len(binBorders);
    const auto sizeTypes = bpy::len(types);
    std::vector<std::string> typesVec {};
    typesVec.reserve((unsigned long) sizeTypes);
    std::vector<double> binBordersVec {};
    binBordersVec.reserve((unsigned long) sizeBorders);
    for(auto i = 0; i < sizeBorders; ++i) {
        binBordersVec.push_back(bpy::extract<double>(binBorders[i]));
    }
    for(auto i = 0; i < sizeTypes; ++i) {
        typesVec.push_back(bpy::extract<std::string>(types[i]));
    }
    auto pyFun = readdy::py::PyFunction<void(readdy::model::HistogramAlongAxisObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::HistogramAlongAxisObservable>(std::move(pyFun), stride, binBordersVec, typesVec, axis);

}

boost::uuids::uuid registerObservable_NParticlesTypes(sim &self, unsigned int stride, bpy::list types, const bpy::object &callbackFun) {
    const auto sizeTypes = bpy::len(types);
    std::vector<std::string> typesVec {};
    typesVec.reserve((unsigned long) sizeTypes);
    for(auto i = 0; i < sizeTypes; ++i) {
        typesVec.push_back(bpy::extract<std::string>(types[i]));
    }
    auto pyFun = readdy::py::PyFunction<void(readdy::model::NParticlesObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::NParticlesObservable>(std::move(pyFun), stride, typesVec);
}

boost::uuids::uuid registerObservable_NParticles(sim &self, unsigned int stride, const bpy::object &callbackFun) {
    auto pyFun = readdy::py::PyFunction<void(readdy::model::NParticlesObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::NParticlesObservable>(std::move(pyFun), stride);
}

#if PY_MAJOR_VERSION >= 3
int
#else
void
#endif
init_numpy()
{
    if(PyArray_API == NULL)
    {
        import_array();
    }
}

// module
BOOST_PYTHON_MODULE (api) {

    init_numpy();
    PyEval_InitThreads();

    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    bpy::docstring_options doc_options;
    doc_options.enable_all();
    bpy::class_<sim, boost::noncopyable>("Simulation")
            .add_property("kbt", &sim::getKBT, &sim::setKBT)
            .add_property("periodic_boundary", &getPeriodicBoundarySimulationWrapper,
                          &setPeriodicBoundarySimulationWrapper)
            .add_property("box_size", &sim::getBoxSize, &setBoxSize)
            .def("register_particle_type", &sim::registerParticleType)
            .def("add_particle", &addParticle)
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
            .def("register_reaction_conversion", &sim::registerConversionReaction, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("register_reaction_enzymatic", &sim::registerEnzymaticReaction, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("register_reaction_fission", &sim::registerFissionReaction, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("register_reaction_fusion", &sim::registerFusionReaction, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("register_reaction_decay", &sim::registerDeathReaction, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("get_recommended_time_step", &sim::getRecommendedTimeStep)
            .def("set_kernel", &sim::setKernel)
            .def("run", &sim::run);

    bpy::class_<kp, boost::noncopyable>("KernelProvider", bpy::no_init)
            .def("get", &kp::getInstance, bpy::return_value_policy<bpy::reference_existing_object>())
            .staticmethod("get")
            .def("load_from_dir", &kp::loadKernelsFromDirectory);

    bpy::class_<vec>("Vec", bpy::init<double, double, double>())
            .def(bpy::self + bpy::self)
            .def(bpy::self - bpy::self)
            .def(double() * bpy::self)
            .def(bpy::self / double())
            .def(bpy::self += bpy::self)
            .def(bpy::self *= double())
            .def(bpy::self == bpy::self)
            .def(bpy::self != bpy::self)
            .def(bpy::self * bpy::self)
            .def(bpy::self_ns::str(bpy::self))
            .def("__getitem__", &pyVec3Bracket);

    bpy::class_<pot2>("Pot2", bpy::init<std::string, boost::python::object, boost::python::object>())
            .def("calc_energy", &pot2::calculateEnergy)
            .def("calc_force", &pot2::calculateForce);

    bpy::class_<kern, boost::noncopyable>("Kernel", bpy::no_init)
            .def("get_name", &kern::getName, bpy::return_value_policy<bpy::reference_existing_object>());

    bpy::class_<uuid>("uuid", bpy::no_init)
            .def("__str__", +[](const uuid& uuid) { return boost::uuids::to_string(uuid);});
}
