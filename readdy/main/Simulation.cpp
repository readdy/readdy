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
    std::vector<std::unique_ptr<readdy::model::potentials::PotentialOrder1>> createdPotentialsOrder1 {};
    std::vector<std::unique_ptr<readdy::model::potentials::PotentialOrder2>> createdPotentialsOrder2 {};
    std::vector<std::unique_ptr<readdy::model::ObservableBase>> createdObservables {};
    std::vector<boost::signals2::connection> observableConnections {};
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

const boost::uuids::uuid& Simulation::registerPotentialOrder2(std::string potential, const std::string &type1, const std::string &type2) {
    ensureKernelSelected();
    // ensure sufficient lifetime
    pimpl->createdPotentialsOrder2.push_back(pimpl->kernel->createPotentialAs<readdy::model::potentials::PotentialOrder2>(potential));
    // add to context
    pimpl->kernel->getKernelContext().registerOrder2Potential(**(pimpl->createdPotentialsOrder2.end() - 1), type1, type2);
    return pimpl->createdPotentialsOrder2.back()->getId();
}

void Simulation::registerPotentialOrder2(readdy::model::potentials::PotentialOrder2 &potential, const std::string &type1, const std::string &type2) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().registerOrder2Potential(potential, type1, type2);
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

void Simulation::registerObservable(const std::string &name, unsigned int stride) {
    ensureKernelSelected();
    pimpl->createdObservables.push_back(pimpl->kernel->createObservable(name));
    (**(pimpl->createdObservables.end()-1)).setStride(stride);
    auto&& connection = pimpl->kernel->registerObservable((*(pimpl->createdObservables.end()-1)).get());
    pimpl->observableConnections.push_back(connection);
}

void Simulation::registerObservable(readdy::model::ObservableBase &observable) {
    ensureKernelSelected();
    auto&& connection = pimpl->kernel->registerObservable(&observable);
    pimpl->observableConnections.push_back(connection);
}


Simulation &Simulation::operator=(Simulation &&rhs) = default;

Simulation::Simulation(Simulation &&rhs) = default;

Simulation::~Simulation() {
    for(auto&& connection : pimpl->observableConnections) {
        connection.disconnect();
    }
}

const boost::uuids::uuid& Simulation::registerPotentialOrder1(std::string name, const std::string &type) {
    ensureKernelSelected();
    pimpl->createdPotentialsOrder1.push_back(pimpl->kernel->createPotentialAs<readdy::model::potentials::PotentialOrder1>(name));
    pimpl->kernel->getKernelContext().registerOrder1Potential(*pimpl->createdPotentialsOrder1.back(), type);
    return pimpl->createdPotentialsOrder1.back()->getId();
}

void Simulation::registerPotentialOrder1(readdy::model::potentials::PotentialOrder1 &potential, const std::string &type) {
    pimpl->kernel->getKernelContext().registerOrder1Potential(potential, type);
}

struct delete_potentialO1 : public std::unary_function<const std::unique_ptr<readdy::model::potentials::PotentialOrder1>, bool> {
    boost::uuids::uuid id;
    delete_potentialO1(const boost::uuids::uuid &id) : id(id) { }
    bool operator()(const std::unique_ptr<readdy::model::potentials::PotentialOrder1> &potential) {
        return id == potential->getId();
    }
};

struct delete_potentialO2 : public std::unary_function<const std::unique_ptr<readdy::model::potentials::PotentialOrder2>, bool> {
    boost::uuids::uuid id;
    delete_potentialO2(const boost::uuids::uuid &id) : id(id) { }
    bool operator()(const std::unique_ptr<readdy::model::potentials::PotentialOrder2> &potential) {
        return id == potential->getId();
    }
};


void Simulation::deregisterPotential(const boost::uuids::uuid &uuid) {
    std::remove_if(pimpl->createdPotentialsOrder1.begin(), pimpl->createdPotentialsOrder1.end(), delete_potentialO1(uuid));
    std::remove_if(pimpl->createdPotentialsOrder2.begin(), pimpl->createdPotentialsOrder2.end(), delete_potentialO2(uuid));
    pimpl->kernel->getKernelContext().deregisterPotential(uuid);
};


NoKernelSelectedException::NoKernelSelectedException(const std::string &__arg) : runtime_error(__arg) { };

