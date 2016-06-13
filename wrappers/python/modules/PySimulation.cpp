#ifdef READDY_WITH_PYTHON

#include <Python.h>
#include <boost/python.hpp>
#include <readdy/Simulation.h>
#include <readdy/plugin/KernelProvider.h>
#include "PotentialWrapper.h"

BOOST_PYTHON_MODULE (simulation) {
    namespace py = boost::python;
    using sim = readdy::Simulation;
    using kp = readdy::plugin::KernelProvider;
    using vec = readdy::model::Vec3;
    using pot2 = readdy::py::PotentialOrder2Wrapper;
    PyEval_InitThreads();
    py::class_<sim, boost::noncopyable>("Simulation")
            .add_property("kbt", &sim::getKBT, &sim::setKBT)
            .add_property("periodic_boundary", &sim::getPeriodicBoundary, &sim::setPeriodicBoundary)
            .add_property("box_size", &sim::getBoxSize, &sim::setBoxSize)
            .def("setKernel", &sim::setKernel);
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