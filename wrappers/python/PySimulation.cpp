#include <readdy/Simulation.h>

#ifdef READDY_WITH_PYTHON
#include <boost/python.hpp>

BOOST_PYTHON_MODULE(simulation) {
        using namespace boost::python;
        class_<readdy::Simulation>("Simulation")
        .def("getKBT", &readdy::Simulation::getKBT)
        .def("setKBT", &readdy::Simulation::setKBT);
}
#endif