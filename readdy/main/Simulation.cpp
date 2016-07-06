//
// Created by Moritz Hoffmann on 18/02/16.
//
#include <readdy/Simulation.h>
#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/programs/Programs.h>
#include <readdy/model/potentials/PotentialsOrder2.h>
#include <readdy/model/potentials/PotentialsOrder1.h>

using namespace readdy;

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
        auto &&diffuseProgram = pimpl->kernel->createProgram<readdy::model::programs::DiffuseProgram>();
        auto &&updateModelProgram = pimpl->kernel->createProgram<readdy::model::programs::UpdateStateModelProgram>();
        auto &&reactionsProgram = pimpl->kernel->createProgram<readdy::model::programs::DefaultReactionProgram>();
        pimpl->kernel->getKernelContext().configure();
        for (readdy::model::time_step_type &&t = 0; t < steps; ++t) {
            updateModelProgram->configure(t, true);
            updateModelProgram->execute();
            diffuseProgram->execute();

            updateModelProgram->configure(t, false);
            updateModelProgram->execute();
            reactionsProgram->execute();
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
    if (fabs(x) <= .5*s[0] && fabs(y) <= .5*s[1] && fabs(z) <= .5*s[2]) {
        readdy::model::Particle p{x, y, z, pimpl->kernel->getKernelContext().getParticleTypeID(type)};
        pimpl->kernel->getKernelStateModel().addParticle(p);
    } else {
        BOOST_LOG_TRIVIAL(error) << "particle position was not in bounds of the simulation box!";
    }

}

void Simulation::registerParticleType(const std::string &name, const double diffusionCoefficient, const double radius) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().setDiffusionConstant(name, diffusionCoefficient);
    pimpl->kernel->getKernelContext().setParticleRadius(name, radius);
}

const std::vector<readdy::model::Vec3> Simulation::getAllParticlePositions() const {
    ensureKernelSelected();
    return pimpl->kernel->getKernelStateModel().getParticlePositions();
}

void Simulation::registerPotentialOrder1(readdy::model::potentials::PotentialOrder1 const* const ptr, const std::string &type) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().registerOrder1Potential(ptr, type);
}

void Simulation::deregisterPotential(const boost::uuids::uuid &uuid) {
    pimpl->kernel->getKernelContext().deregisterPotential(uuid);
};

boost::uuids::uuid Simulation::registerHarmonicRepulsionPotential(std::string particleTypeA, std::string particleTypeB, double forceConstant) {
    ensureKernelSelected();
    auto ptr = pimpl->kernel->createPotentialAs<readdy::model::potentials::HarmonicRepulsion>();
    ptr->setForceConstant(forceConstant);
    return pimpl->kernel->getKernelContext().registerOrder2Potential(ptr.get(), particleTypeA, particleTypeB);
}

boost::uuids::uuid Simulation::registerWeakInteractionPiecewiseHarmonicPotential(std::string particleTypeA, std::string particleTypeB, double forceConstant, double desiredParticleDistance,
                                                                                 double depth, double noInteractionDistance) {
    ensureKernelSelected();
    auto ptr = pimpl->kernel->createPotentialAs<readdy::model::potentials::WeakInteractionPiecewiseHarmonic>();
    ptr->setForceConstant(forceConstant);
    ptr->setDesiredParticleDistance(desiredParticleDistance);
    ptr->setDepthAtDesiredDistance(depth);
    ptr->setNoInteractionDistance(noInteractionDistance);
    return pimpl->kernel->getKernelContext().registerOrder2Potential(ptr.get(), particleTypeA, particleTypeB);
}

boost::uuids::uuid Simulation::registerBoxPotential(std::string particleType, double forceConstant, readdy::model::Vec3 origin, readdy::model::Vec3 extent, bool considerParticleRadius) {
    ensureKernelSelected();
    auto ptr = pimpl->kernel->createPotentialAs<readdy::model::potentials::CubePotential>();
    ptr->setOrigin(origin);
    ptr->setExtent(extent);
    ptr->setConsiderParticleRadius(considerParticleRadius);
    ptr->setForceConstant(forceConstant);
    return pimpl->kernel->getKernelContext().registerOrder1Potential(ptr.get(), particleType);
}

void Simulation::registerPotentialOrder2(model::potentials::PotentialOrder2 const* const ptr, const std::string &type1, const std::string &type2) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().registerOrder2Potential(ptr, type1, type2);
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

boost::uuids::uuid Simulation::registerObservable(readdy::model::ObservableBase &observable) {
    ensureKernelSelected();
    boost::uuids::random_generator uuid_gen;
    auto uuid = uuid_gen();
    auto&& connection = pimpl->kernel->connectObservable(&observable);
    pimpl->observableConnections.emplace(uuid, std::move(connection));
    return uuid;
}

void Simulation::deregisterObservable(const boost::uuids::uuid uuid) {
    pimpl->observableConnections.erase(uuid);
    if(pimpl->observables.find(uuid) != pimpl->observables.end()) {
        pimpl->observables.erase(uuid);
    }
}

std::vector<std::string> Simulation::getAvailableObservables() {
    ensureKernelSelected();
    // TODO compile a list of observables
    return {"hallo"};
}


Simulation &Simulation::operator=(Simulation &&rhs) = default;

Simulation::Simulation(Simulation &&rhs) = default;

Simulation::~Simulation() = default;

const boost::uuids::uuid &Simulation::registerConversionReaction(const std::string &name, const std::string &from, const std::string &to, const double &rate) {
    ensureKernelSelected();
    return pimpl->kernel->getKernelContext().registerConversionReaction(name, from, to, rate);
}

const boost::uuids::uuid &Simulation::registerEnzymaticReaction(const std::string &name, const std::string &catalyst, const std::string &from, const std::string &to, const double &rate,
                                                                const double &eductDistance) {
    ensureKernelSelected();
    return pimpl->kernel->getKernelContext().registerEnzymaticReaction(name, catalyst, from, to, rate, eductDistance);
}

const boost::uuids::uuid &Simulation::registerFissionReaction(const std::string &name, const std::string &from, const std::string &to1, const std::string &to2, const double productDistance,
                                                              const double &rate, const double &weight1, const double &weight2) {
    ensureKernelSelected();
    return pimpl->kernel->getKernelContext().registerFissionReaction(name, from, to1, to2, productDistance, rate, weight1, weight2);
}

const boost::uuids::uuid &Simulation::registerFusionReaction(const std::string &name, const std::string &from1, const std::string &from2, const std::string &to, const double &rate,
                                                             const double &eductDistance, const double &weight1, const double &weight2) {
    ensureKernelSelected();
    return pimpl->kernel->getKernelContext().registerFusionReaction(name, from1, from2, to, rate, eductDistance, weight1, weight2);
}

const boost::uuids::uuid &Simulation::registerDeathReaction(const std::string &name, const std::string &particleType, const double &rate) {
    ensureKernelSelected();
    return pimpl->kernel->getKernelContext().registerDeathReaction(name, particleType, rate);
}

std::vector<readdy::model::Vec3> Simulation::getParticlePositions(std::string type) {
    unsigned int typeId = pimpl->kernel->getKernelContext().getParticleTypeID(type);
    const auto particles = pimpl->kernel->getKernelStateModel().getParticles();
    std::vector<readdy::model::Vec3> positions;
    for(auto&& p : particles) {
        if(p.getType() == typeId) {
            positions.push_back(p.getPos());
        }
    }
    return positions;
}

/*template boost::uuids::uuid Simulation::registerObservable<readdy::model::ParticlePositionObservable>(std::function<void(const typename readdy::model::ParticlePositionObservable::result_t)>, unsigned int);
template boost::uuids::uuid Simulation::registerObservable<readdy::model::RadialDistributionObservable>(
        std::function<void(const typename readdy::model::RadialDistributionObservable::result_t)>, unsigned int,
        std::vector<double>, const std::string, const std::string, double
);
template boost::uuids::uuid Simulation::registerObservable<readdy::model::CenterOfMassObservable>(
        std::function<void(const typename readdy::model::CenterOfMassObservable::result_t)>, unsigned int, const std::vector<std::string>
);
template boost::uuids::uuid Simulation::registerObservable<readdy::model::HistogramAlongAxisObservable>(
        std::function<void(const typename readdy::model::HistogramAlongAxisObservable::result_t)>, unsigned int, std::vector<double> , std::vector<std::string> , unsigned int
);
template boost::uuids::uuid Simulation::registerObservable<readdy::model::NParticlesObservable>(
        std::function<void(const typename readdy::model::NParticlesObservable::result_t)>, unsigned int
);*/


NoKernelSelectedException::NoKernelSelectedException(const std::string &__arg) : runtime_error(__arg) { };

