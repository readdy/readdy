#include <readdy/Simulation.h>

#ifdef READDY_WITH_PYTHON

#include <Python.h>
#include <boost/python.hpp>

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#pragma message "content of PYTHON_HEX: " STR(PY_VERSION_HEX)

BOOST_PYTHON_MODULE(simulation) {
        using namespace boost::python;
        class_<readdy::Simulation>("Simulation")
        .def("getKBT", &readdy::Simulation::getKBT)
        .def("setKBT", &readdy::Simulation::setKBT);
}
#endif