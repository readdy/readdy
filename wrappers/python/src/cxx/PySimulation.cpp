
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include "PyConverters.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <readdy/Simulation.h>
#include <readdy/plugin/KernelProvider.h>
#include "PyPotential.h"
#include "PyFunction.h"

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
    std::vector<std::string> typesVec {};
    const auto len = boost::python::len(types);
    typesVec.reserve((unsigned long) len);
    for(auto i = 0; i < len; ++i) {
        typesVec.push_back(boost::python::extract<std::string>(types[i]));
    }
    auto pyFun = readdy::py::PyFunction<void(readdy::model::CenterOfMassObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::CenterOfMassObservable>(std::move(pyFun), stride, typesVec);
}

#if PY_MAJOR_VERSION >= 3
int
#else
void
#endif
init_numpy()
{
    import_array();
}

// module
BOOST_PYTHON_MODULE (simulation) {

    PyEval_InitThreads();

    init_numpy();

    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    readdy::py::std_vector_to_python_converter<double>();
    readdy::py::std_pair_to_python_converter<std::vector<double>, std::vector<double>>();

    bpy::class_<sim, boost::noncopyable>("Simulation")
            .add_property("kbt", &sim::getKBT, &sim::setKBT)
            .add_property("periodic_boundary", &getPeriodicBoundarySimulationWrapper, &setPeriodicBoundarySimulationWrapper)
            .add_property("box_size", &sim::getBoxSize, &setBoxSize)
            .def("registerParticleType", &sim::registerParticleType)
            .def("addParticle", &addParticle)
            .def("isKernelSelected", &sim::isKernelSelected)
            .def("getSelectedKernelType", &getSelectedKernelType)
            .def("registerPotentialOrder2", &registerPotentialOrder2)
            .def("registerHarmonicRepulsionPotential", &sim::registerHarmonicRepulsionPotential)
            .def("registerWeakInteractionPiecewiseHarmonicPotential", &sim::registerWeakInteractionPiecewiseHarmonicPotential)
            .def("registerBoxPotential", &sim::registerBoxPotential)
            .def("getParticlePositions", &sim::getParticlePositions)
            .def("registerObservable_ParticlePositions", &registerObservable_ParticlePositions)
            .def("registerObservable_RadialDistribution", &registerObservable_RadialDistribution)
            .def("registerObservable_CenterOfMass", &registerObservable_CenterOfMass)
            .def("registerConversionReaction", &sim::registerConversionReaction, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("registerEnzymaticReaction", &sim::registerEnzymaticReaction, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("registerFissionReaction", &sim::registerFissionReaction, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("registerFusionReaction", &sim::registerFusionReaction, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("registerDeathReaction", &sim::registerDeathReaction, bpy::return_value_policy<bpy::reference_existing_object>())
            .def("setKernel", &sim::setKernel)
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

    bpy::class_<std::vector<vec>>("Vecvec")
            .def(boost::python::vector_indexing_suite<std::vector<vec>>());

    bpy::class_<pot2>("Pot2", bpy::init<std::string, boost::python::object, boost::python::object>())
            .def("calc_energy", &pot2::calculateEnergy)
            .def("calc_force", &pot2::calculateForce);

    bpy::class_<kern, boost::noncopyable>("Kernel", bpy::no_init)
            .def("getName", &kern::getName, bpy::return_value_policy<bpy::reference_existing_object>());

    bpy::class_<uuid>("uuid", bpy::no_init)
            .def("__str__", +[](const uuid& uuid) { return boost::uuids::to_string(uuid);});
}
