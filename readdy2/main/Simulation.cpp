//
// Created by Moritz Hoffmann on 18/02/16.
//
#include <readdy/Simulation.h>
#include <readdy/common/make_unique.h>
#include <readdy/plugin/Kernel.h>
#include <readdy/plugin/KernelProvider.h>

using namespace readdy;

struct Simulation::Impl {
    std::shared_ptr<readdy::plugin::Kernel> kernel = nullptr;
};

double Simulation::getKBT() const {
    if (isKernelSelected()) {
        return pimpl->kernel->getKernelContext().getKBT();
    }
    throw NoKernelSelectedException("No Kernel was set");
}

void Simulation::setKBT(double kBT) {
    if (isKernelSelected()) {
        pimpl->kernel->getKernelContext().setKBT(kBT);
    } else {
        throw NoKernelSelectedException("No Kernel was set");
    }
}

void Simulation::setBoxSize(double dx, double dy, double dz) {
    if (isKernelSelected()) {
        pimpl->kernel->getKernelContext().setBoxSize(dx, dy, dz);
    } else {
        throw NoKernelSelectedException("no kernel was set");
    }
}

void Simulation::setPeriodicBoundary(std::array<bool, 3> periodic) {
    if(isKernelSelected()) {
        pimpl->kernel->getKernelContext().setPeriodicBoundary(periodic[0], periodic[1], periodic[2]);
    } else {
        throw NoKernelSelectedException("no kernel was set");
    }
}

Simulation::Simulation() : pimpl(std::make_unique<Simulation::Impl>()) { }

std::array<double, 3> Simulation::getBoxSize() const {
    if (isKernelSelected()) {
        return pimpl->kernel->getKernelContext().getBoxSize();
    }
    throw NoKernelSelectedException("No Kernel was set");
}

std::array<bool, 3> Simulation::getPeriodicBoundary() const {
    if (isKernelSelected()) {
        return pimpl->kernel->getKernelContext().getPeriodicBoundary();
    }
    throw NoKernelSelectedException("No Kernel was set");
}

Simulation::Simulation(const Simulation &rhs) : pimpl(std::make_unique<Simulation::Impl>(*rhs.pimpl)) { }

Simulation &Simulation::operator=(const Simulation &rhs) {
    *pimpl = *rhs.pimpl;
    return *this;
}


void Simulation::run(const readdy::model::time_step_type steps, const double timeStep) {
    if (pimpl->kernel) {
        const auto& kernel = pimpl->kernel;
        {
            BOOST_LOG_TRIVIAL(debug) << "available programs: ";
            for (auto &&p : pimpl->kernel->getAvailablePrograms()) {
                BOOST_LOG_TRIVIAL(debug) << "\t" << p;
            }
        }
        kernel->getKernelContext().setTimeStep(timeStep);
        {
            auto&& diffuseProgram = kernel->createProgram("Diffuse");
            for (auto&& t = 0; t < steps; ++t) {
                diffuseProgram->execute();
            }
        }
    } else {
        BOOST_LOG_TRIVIAL(error) << "no kernel was selected!";
        throw NoKernelSelectedException("Before running the simulation, you must select a kernel");
    }
}

void Simulation::setKernel(const std::string& kernel) {
    if (pimpl->kernel) {
        BOOST_LOG_TRIVIAL(debug) << "replacing kernel \"" << pimpl->kernel->getName() << "\" with \"" << kernel << "\"";
    }
    pimpl->kernel = readdy::plugin::KernelProvider::getInstance().get(kernel);
}

bool Simulation::isKernelSelected() const {
    return pimpl->kernel ? true : false;
}

const std::string& Simulation::getSelectedKernelType() const {
    if (isKernelSelected()) {
        return pimpl->kernel->getName();
    }
    throw NoKernelSelectedException("You must select a kernel by setKernel befor asking for its type");
}

void Simulation::addParticle(double x, double y, double z, const std::string& type) {
    if (isKernelSelected()) {
        const auto&& s = getBoxSize();
        if (0 <= x && x <= s[0] && 0 <= y && y <= s[1] && 0 <= z && z <= s[2]) {
            readdy::model::Particle p{x, y, z, pimpl->kernel->getKernelContext().getParticleTypeID(type)};
            pimpl->kernel->getKernelStateModel().addParticle(p);
        } else {
            BOOST_LOG_TRIVIAL(error) << "particle position was not in bounds of the simulation box!";
        }
    } else {
        throw NoKernelSelectedException("no kernel was set!");
    }
}

void Simulation::registerParticleType(const std::string name, const double diffusionCoefficient) {
    if(isKernelSelected()) {
        pimpl->kernel->getKernelContext().setDiffusionConstant(name, diffusionCoefficient);
    } else {
        throw NoKernelSelectedException("no kernel was selected");
    }
}

const std::vector<readdy::model::Vec3> Simulation::getParticlePositions() const{
    return pimpl->kernel->getKernelStateModel().getParticlePositions();
}




Simulation &Simulation::operator=(Simulation &&rhs) = default;

Simulation::Simulation(Simulation &&rhs) = default;

Simulation::~Simulation() = default;


NoKernelSelectedException::NoKernelSelectedException(const std::string &__arg) : runtime_error(__arg) { };

