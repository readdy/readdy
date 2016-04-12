//
// Created by Moritz Hoffmann on 18/02/16.
//
#include <readdy/Simulation.h>
#include <boost/make_unique.hpp>

#include <boost/predef.h>
#ifdef BOOST_OS_MACOS
#include <array>
#endif

using namespace readdy;

struct Simulation::Impl {
    double kBT = 0;
    std::array<double, 3> box_size {};
    std::array<bool, 3> periodic_boundary {};
};

double Simulation::getKBT() const {
    return (*pimpl).kBT;
}

void Simulation::setKBT(double kBT) {
    (*pimpl).kBT = kBT;
}

void Simulation::setBoxSize(double dx, double dy, double dz) {
    (*pimpl).box_size = {dx, dy, dz};
}

void Simulation::setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z) {
    (*pimpl).periodic_boundary = {pb_x, pb_y, pb_z};
}

Simulation::Simulation() : pimpl(boost::make_unique<Simulation::Impl>()) {}

std::array<double, 3> Simulation::getBoxSize() const {
    return std::array<double, 3>((*pimpl).box_size);
}

std::array<bool, 3> Simulation::getPeriodicBoundary() const {
    return std::array<bool, 3>((*pimpl).periodic_boundary);
}

Simulation::Simulation(const Simulation &rhs) : pimpl(boost::make_unique<Simulation::Impl>(*rhs.pimpl)){}

Simulation &Simulation::operator=(const Simulation &rhs) {
    *pimpl = *rhs.pimpl;
    return *this;
}


Simulation& Simulation::operator=(Simulation &&rhs) = default;
Simulation::Simulation(Simulation &&rhs) = default;
Simulation::~Simulation() = default;



