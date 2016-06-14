#ifdef READDY_WITH_PYTHON

#include <Python.h>
#include <boost/python.hpp>
#include <readdy/Simulation.h>
#include <readdy/plugin/KernelProvider.h>
#include "PyPotential.h"

namespace py = boost::python;
using sim = readdy::Simulation;
using kp = readdy::plugin::KernelProvider;
using vec = readdy::model::Vec3;
using pot2 = readdy::py::PotentialOrder2Wrapper;
using model = readdy::model::KernelStateModel;
using ctx = readdy::model::KernelContext;

boost::python::list getPeriodicBoundarySimulationWrapper(const sim &simulation) {
    auto p = simulation.getPeriodicBoundary();
    py::list result;
    result.append<bool>(p[0]);
    result.append<bool>(p[1]);
    result.append<bool>(p[2]);
    return result;
}

void setPeriodicBoundarySimulationWrapper(sim &self, py::list list) {
    self.setPeriodicBoundary({py::extract<bool>(list[0]), py::extract<bool>(list[1]), py::extract<bool>(list[2])});
}

// thin wrappers
void setBoxSize(sim &self, const vec& size) { /* explicitely choose void(vec) signature */ self.setBoxSize(size); }
std::string getSelectedKernelType(sim &self) { /* discard const reference */ return self.getSelectedKernelType(); }
void addParticle(sim& self, const std::string& type, const vec& pos) { self.addParticle(pos[0], pos[1], pos[2], type); }
void registerPotentialOrder2(sim& self, pot2& potential, std::string type1, std::string type2) { self.registerPotentialOrder2(potential, type1, type2); }
void registerPotentialOrder2_name(sim& self, std::string potentialType, std::string type1, std::string type2) { self.registerPotentialOrder2(potentialType, type1, type2); }

// module
BOOST_PYTHON_MODULE (simulation) {

    PyEval_InitThreads();
    py::class_<sim, boost::noncopyable>("Simulation")
            .add_property("kbt", &sim::getKBT, &sim::setKBT)
            .add_property("periodic_boundary", &getPeriodicBoundarySimulationWrapper, &setPeriodicBoundarySimulationWrapper)
            .add_property("box_size", &sim::getBoxSize, &setBoxSize)
            .def("registerParticleType", &sim::registerParticleType)
            .def("addParticle", &addParticle)
            .def("isKernelSelected", &sim::isKernelSelected)
            .def("getSelectedKernelType", &getSelectedKernelType)
            .def("registerPotentialOrder2", &registerPotentialOrder2)
            .def("registerPotentialOrder2", &registerPotentialOrder2_name)
            .def("setKernel", &sim::setKernel)
            .def("run", &sim::run);
    py::class_<kp, boost::noncopyable>("KernelProvider", py::no_init)
            .def("get", &kp::getInstance, py::return_value_policy<py::reference_existing_object>())
            .staticmethod("get")
            .def("load_from_dir", &kp::loadKernelsFromDirectory);
    py::class_<vec>("Vec", py::init<double, double, double>())
            .def(py::self + py::self)
            .def(py::self - py::self)
            .def(double() * py::self)
            .def(py::self += py::self)
            .def(py::self *= double())
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def(py::self * py::self)
            .def(py::self_ns::str(py::self))
            .def("__getitem__", &vec::operator[]);

    py::class_<pot2>("Pot2", py::init<std::string, boost::python::object, boost::python::object>())
            .def("calc_energy", &pot2::calculateEnergy)
            .def("calc_force", &pot2::calculateForce);
}

#endif