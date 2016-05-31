#ifdef READDY_WITH_PYTHON

#include <Python.h>
#include <boost/python.hpp>
#include <readdy/Simulation.h>
#include <readdy/plugin/KernelProvider.h>

BOOST_PYTHON_MODULE (simulation) {
    namespace py = boost::python;
    using sim = readdy::Simulation;
    using kp = readdy::plugin::KernelProvider;
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
}

#endif