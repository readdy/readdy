#include <readdy/Simulation.h>

#ifdef READDY_WITH_PYTHON

#include <Python.h>
#include <boost/python.hpp>

BOOST_PYTHON_MODULE (simulation) {
    namespace py = boost::python;
    using sim = readdy::Simulation;
    PyEval_InitThreads();
    py::class_<sim, boost::noncopyable>("Simulation")
            .add_property("kbt", &sim::getKBT, &sim::setKBT)
            .add_property("periodic_boundary", &sim::getPeriodicBoundary, &sim::setPeriodicBoundary)
            .add_property("box_size", &sim::getBoxSize, &sim::setBoxSize);
}

#endif