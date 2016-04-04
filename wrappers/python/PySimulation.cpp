#include <readdy/Simulation.h>

#ifdef READDY_WITH_PYTHON

#include <Python.h>
#include <boost/python.hpp>

#if PY_MAJOR_VERSION >= 3
#if PY_VERSION_HEX >= 0x03000000
#warning "definitely using python 3 !!! !! ! !"
#else
#error "definitely something wrong :-("
#endif
#endif

BOOST_PYTHON_MODULE(simulation) {
        using namespace boost::python;
        class_<readdy::Simulation>("Simulation")
        .def("getKBT", &readdy::Simulation::getKBT)
        .def("setKBT", &readdy::Simulation::setKBT);
}
#endif