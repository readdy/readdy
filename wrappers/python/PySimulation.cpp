#include <Simulation.h>

#ifdef READDY_WITH_PYTHON
#include <boost/python.hpp>

BOOST_PYTHON_MODULE(simulation) {
        using namespace boost::python;
        class_<ReaDDy::Simulation>("Simulation")
        .def("getKBT", &ReaDDy::Simulation::getKBT)
        .def("setKBT", &ReaDDy::Simulation::setKBT);
}
#endif