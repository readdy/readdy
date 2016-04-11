#include <readdy/Simulation.h>

#ifdef READDY_WITH_PYTHON

#include <Python.h>
#include <boost/python.hpp>

BOOST_PYTHON_MODULE (simulation) {
    namespace py = boost::python;
    using sim = readdy::Simulation;
    py::class_<sim>("Simulation")
            .def("getKBT", &sim::getKBT)
            .def("setKBT", &sim::setKBT)
            .def("setPeriodicBoundary", &sim::setPeriodicBoundary)
            .def("getPeriodicBoundary", &sim::getPeriodicBoundary)
            .def("setBoxSize", &sim::setBoxSize)
            .def("getBoxSize", &sim::getBoxSize);
}

#endif