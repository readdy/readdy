//
// Created by Moritz Hoffmann on 18/02/16.
//
#include <readdy/Simulation.h>
#include <boost/make_unique.hpp>
#include <readdy/plugin/Kernel.h>

using namespace readdy;

struct Simulation::Impl {
    std::shared_ptr<readdy::plugin::Kernel> kernel = nullptr;
};

double Simulation::getKBT() const {
    if(isKernelSet()) {
        return pimpl->kernel->getKernelContext()->getKBT();
    }
    throw NoKernelSelectedException("No Kernel was set");
}

void Simulation::setKBT(double kBT) {
    if(isKernelSet()) {
        pimpl->kernel->getKernelContext()->setKBT(kBT);
    } else {
        throw NoKernelSelectedException("No Kernel was set");
    }
}

void Simulation::setBoxSize(double dx, double dy, double dz) {
    if(isKernelSet()) {
        pimpl->kernel->getKernelContext()->setBoxSize(dx, dy, dz);
    } else {
        throw NoKernelSelectedException("no kernel was set");
    }
}

void Simulation::setPeriodicBoundary(bool pb_x, bool pb_y, bool pb_z) {
    if(isKernelSet()) {
        pimpl->kernel->getKernelContext()->setPeriodicBoundary(pb_x, pb_y, pb_z);
    } else {
        throw NoKernelSelectedException("no kernel was set");
    }
}

Simulation::Simulation() : pimpl(boost::make_unique<Simulation::Impl>()) {}

std::array<double, 3> Simulation::getBoxSize() const {
    if(isKernelSet()) {
        return pimpl->kernel->getKernelContext()->getBoxSize();
    }
    throw NoKernelSelectedException("No Kernel was set");
}

std::array<bool, 3> Simulation::getPeriodicBoundary() const {
    if(isKernelSet()) {
        return pimpl->kernel->getKernelContext()->getPeriodicBoundary();
    }
    throw NoKernelSelectedException("No Kernel was set");
}

Simulation::Simulation(const Simulation &rhs) : pimpl(boost::make_unique<Simulation::Impl>(*rhs.pimpl)){}

Simulation &Simulation::operator=(const Simulation &rhs) {
    *pimpl = *rhs.pimpl;
    return *this;
}



void Simulation::run(const size_t steps, const double timeStep) {
    if(pimpl->kernel) {
        auto kernel = pimpl->kernel;
        {
            BOOST_LOG_TRIVIAL(debug) << "available programs: ";
            for (auto &&p : pimpl->kernel->getAvailablePrograms()) {
                BOOST_LOG_TRIVIAL(debug) << "\t" << p;
            }
        }
        kernel->getKernelContext()->setTimeStep(timeStep);
        for (int t = 0; t < steps; ++t) {

        }
    } else {
        BOOST_LOG_TRIVIAL(error) << "no kernel was selected!";
        throw NoKernelSelectedException("Before running the simulation, you must select a kernel");
    }
}

void Simulation::setKernel(const std::string kernel) {
    if(pimpl->kernel) {
        BOOST_LOG_TRIVIAL(debug) << "replacing kernel \""<< pimpl->kernel->getName() << "\" with \"" << kernel << "\"";
    }
    pimpl->kernel = readdy::plugin::KernelProvider::getInstance().get(kernel);
}

bool Simulation::isKernelSet() const {
    return pimpl->kernel ? true : false;
}

std::string Simulation::getSelectedKernelType() const {
    if(isKernelSet()) {
        return pimpl->kernel->getName();
    }
    throw NoKernelSelectedException("You must select a kernel by setKernel befor asking for its type");
}


Simulation& Simulation::operator=(Simulation &&rhs) = default;
Simulation::Simulation(Simulation &&rhs) = default;
Simulation::~Simulation() = default;


NoKernelSelectedException::NoKernelSelectedException(const std::string &__arg) : runtime_error(__arg) {};

