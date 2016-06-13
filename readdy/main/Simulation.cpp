//
// Created by Moritz Hoffmann on 18/02/16.
//
#include <readdy/Simulation.h>
#include <readdy/common/make_unique.h>
#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>

using namespace readdy;

struct Simulation::Impl {
    std::unique_ptr<readdy::model::Kernel> kernel;
    std::vector<std::unique_ptr<readdy::model::potentials::Potential>> createdPotentials {};
};

double Simulation::getKBT() const {
    ensureKernelSelected();
    return pimpl->kernel->getKernelContext().getKBT();

}

void Simulation::setKBT(double kBT) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().setKBT(kBT);

}

void Simulation::setBoxSize(const readdy::model::Vec3 &boxSize) {
    setBoxSize(boxSize[0], boxSize[1], boxSize[2]);
}

void Simulation::setPeriodicBoundary(std::array<bool, 3> periodic) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().setPeriodicBoundary(periodic[0], periodic[1], periodic[2]);

}

Simulation::Simulation() : pimpl(std::make_unique<Simulation::Impl>()) { }

readdy::model::Vec3 Simulation::getBoxSize() const {
    ensureKernelSelected();
    return readdy::model::Vec3(pimpl->kernel->getKernelContext().getBoxSize());

}

std::array<bool, 3> Simulation::getPeriodicBoundary() const {
    ensureKernelSelected();
    return pimpl->kernel->getKernelContext().getPeriodicBoundary();
}

void Simulation::run(const readdy::model::time_step_type steps, const double timeStep) {
    ensureKernelSelected();
    {
        BOOST_LOG_TRIVIAL(debug) << "available programs: ";
        for (auto &&p : pimpl->kernel->getAvailablePrograms()) {
            BOOST_LOG_TRIVIAL(debug) << "\t" << p;
        }
    }
    pimpl->kernel->getKernelContext().setTimeStep(timeStep);
    {
        auto &&diffuseProgram = pimpl->kernel->createProgram("Diffuse");
        for (auto &&t = 0; t < steps; ++t) {
            diffuseProgram->execute();
        }
    }
}

void Simulation::setKernel(const std::string &kernel) {
    if (isKernelSelected()) {
        BOOST_LOG_TRIVIAL(debug) << "replacing kernel \"" << pimpl->kernel->getName() << "\" with \"" << kernel << "\"";
    }
    pimpl->kernel = readdy::plugin::KernelProvider::getInstance().create(kernel);
}

bool Simulation::isKernelSelected() const {
    return pimpl->kernel ? true : false;
}

const std::string &Simulation::getSelectedKernelType() const {
    ensureKernelSelected();
    return pimpl->kernel->getName();
}

void Simulation::addParticle(double x, double y, double z, const std::string &type) {
    ensureKernelSelected();
    const auto &&s = getBoxSize();
    if (0 <= x && x <= s[0] && 0 <= y && y <= s[1] && 0 <= z && z <= s[2]) {
        readdy::model::Particle p{x, y, z, pimpl->kernel->getKernelContext().getParticleTypeID(type)};
        pimpl->kernel->getKernelStateModel().addParticle(p);
    } else {
        BOOST_LOG_TRIVIAL(error) << "particle position was not in bounds of the simulation box!";
    }

}

void Simulation::registerParticleType(const std::string &name, const double diffusionCoefficient) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().setDiffusionConstant(name, diffusionCoefficient);
}

const std::vector<readdy::model::Vec3> Simulation::getParticlePositions() const {
    ensureKernelSelected();
    return pimpl->kernel->getKernelStateModel().getParticlePositions();
}

void Simulation::registerPotentialOrder2(std::string potential, const std::string &type1, const std::string &type2) {
    ensureKernelSelected();
    // ensure sufficient lifetime
    pimpl->createdPotentials.push_back(pimpl->kernel->createPotential(potential));
    // add to context
    pimpl->kernel->getKernelContext().registerPotentialForTypes(**(pimpl->createdPotentials.end()-1), type1, type2);
}

void Simulation::registerPotentialOrder2(readdy::model::potentials::Potential &potential, const std::string &type1, const std::string &type2) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().registerPotentialForTypes(potential, type1, type2);
}

void Simulation::ensureKernelSelected() const {
    if (!isKernelSelected()) {
        throw NoKernelSelectedException("No kernel was selected!");
    }
}

void Simulation::setBoxSize(double dx, double dy, double dz) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().setBoxSize(dx, dy, dz);
}


Simulation &Simulation::operator=(Simulation &&rhs) = default;

Simulation::Simulation(Simulation &&rhs) = default;

Simulation::~Simulation() = default;


NoKernelSelectedException::NoKernelSelectedException(const std::string &__arg) : runtime_error(__arg) { };

