/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


//
// Created by Moritz Hoffmann on 18/02/16.
//
#include <readdy/api/Simulation.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/Utils.h>
#include <readdy/model/observables/io/Trajectory.h>
#include <readdy/io/File.h>

namespace rmr = readdy::model::reactions;
namespace rmp = readdy::model::actions;

namespace readdy {
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

Simulation::Simulation() : pimpl(std::make_unique<Simulation::Impl>()) {}

readdy::model::Vec3 Simulation::getBoxSize() const {
    ensureKernelSelected();
    return readdy::model::Vec3(pimpl->kernel->getKernelContext().getBoxSize());

}

std::array<bool, 3> Simulation::getPeriodicBoundary() const {
    ensureKernelSelected();
    return pimpl->kernel->getKernelContext().getPeriodicBoundary();
}

void Simulation::run(const readdy::model::observables::time_step_type steps, const double timeStep) {
    ensureKernelSelected();
    {
        log::console()->debug("available actions: ");
        for (auto &&p : pimpl->kernel->getAvailableActions()) {
            log::console()->debug("\t {}", p);
        }
    }
    runScheme().configure(timeStep)->run(steps);
}

void Simulation::setKernel(const std::string &kernel) {
    if (isKernelSelected()) {
        log::console()->debug("replacing kernel \"{}\" with \"{}\"", pimpl->kernel->getName(), kernel);
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
    if (fabs(x) <= .5 * s[0] && fabs(y) <= .5 * s[1] && fabs(z) <= .5 * s[2]) {
        readdy::model::Particle p{x, y, z, pimpl->kernel->getKernelContext().getParticleTypeID(type)};
        pimpl->kernel->getKernelStateModel().addParticle(p);
    } else {
        log::console()->error("particle position was not in bounds of the simulation box!");
    }

}

Simulation::particle_t::type_type
Simulation::registerParticleType(const std::string &name, const double diffusionCoefficient, const double radius) {
    ensureKernelSelected();
    auto& context = pimpl->kernel->getKernelContext();
    context.setDiffusionConstant(name, diffusionCoefficient);
    context.setParticleRadius(name, radius);
    return context.getParticleTypeID(name);
}

const std::vector<readdy::model::Vec3> Simulation::getAllParticlePositions() const {
    ensureKernelSelected();
    return pimpl->kernel->getKernelStateModel().getParticlePositions();
}

void Simulation::deregisterPotential(const short uuid) {
    pimpl->kernel->getKernelContext().deregisterPotential(uuid);
};

const short
Simulation::registerHarmonicRepulsionPotential(const std::string& particleTypeA, const std::string& particleTypeB,
                                               double forceConstant) {
    using potential_t = readdy::model::potentials::HarmonicRepulsion;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleTypeA, particleTypeB, forceConstant);
}

const short
Simulation::registerWeakInteractionPiecewiseHarmonicPotential(const std::string& particleTypeA,
                                                              const std::string& particleTypeB, double forceConstant,
                                                              double desiredParticleDistance, double depth,
                                                              double noInteractionDistance) {
    using potential_t = readdy::model::potentials::WeakInteractionPiecewiseHarmonic;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleTypeA, particleTypeB, forceConstant,
                                                         desiredParticleDistance, depth, noInteractionDistance);
}

const short
Simulation::registerBoxPotential(const std::string &particleType, double forceConstant,
                                 const readdy::model::Vec3 &origin, const readdy::model::Vec3 &extent,
                                 bool considerParticleRadius) {
    using potential_t = readdy::model::potentials::CubePotential;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleType, forceConstant, origin, extent,
                                                         considerParticleRadius);
}

const short
Simulation::registerSpherePotential(std::string particleType, double forceConstant, const readdy::model::Vec3 &origin,
                                    double radius) {
    using potential_t = readdy::model::potentials::SpherePotential;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleType, forceConstant, origin, radius);
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

ObservableHandle Simulation::registerObservable(readdy::model::observables::ObservableBase &observable) {
    ensureKernelSelected();
    auto uuid = pimpl->counter++;
    auto &&connection = pimpl->kernel->connectObservable(&observable);
    pimpl->observableConnections.emplace(uuid, std::move(connection));
    return {uuid, nullptr};
}

void Simulation::deregisterObservable(const unsigned long uuid) {
    pimpl->observableConnections.erase(uuid);
    if (pimpl->observables.find(uuid) != pimpl->observables.end()) {
        pimpl->observables.erase(uuid);
    }
}

void Simulation::deregisterObservable(const ObservableHandle &uuid) {
    deregisterObservable(uuid.getId());
}


std::vector<std::string> Simulation::getAvailableObservables() {
    ensureKernelSelected();
    // TODO compile a list of observables
    return {"hallo"};
}


Simulation &Simulation::operator=(Simulation &&rhs) = default;

Simulation::Simulation(Simulation &&rhs) = default;

Simulation::~Simulation() {
    // close trajectory file if present
    closeTrajectoryFile();
};

const short
Simulation::registerConversionReaction(const std::string &name, const std::string &from, const std::string &to,
                                       const double rate) {
    ensureKernelSelected();
    namespace rmr = readdy::model::reactions;
    auto reaction = pimpl->kernel->createConversionReaction(name, from, to, rate);
    return pimpl->kernel->getKernelContext().registerReaction(std::move(reaction));
}

const short
Simulation::registerEnzymaticReaction(const std::string &name, const std::string &catalyst, const std::string &from,
                                      const std::string &to, const double rate,
                                      const double eductDistance) {
    ensureKernelSelected();
    namespace rmr = readdy::model::reactions;
    auto reaction = pimpl->kernel->createEnzymaticReaction(name, catalyst, from, to, rate, eductDistance);
    return pimpl->kernel->getKernelContext().registerReaction(std::move(reaction));
}

const short
Simulation::registerFissionReaction(const std::string &name, const std::string &from, const std::string &to1,
                                    const std::string &to2,
                                    const double rate, const double productDistance, const double weight1,
                                    const double weight2) {
    ensureKernelSelected();
    auto reaction = pimpl->kernel->createFissionReaction(name, from, to1, to2, rate, productDistance, weight1, weight2);
    return pimpl->kernel->getKernelContext().registerReaction(std::move(reaction));
}

const short
Simulation::registerFusionReaction(const std::string &name, const std::string &from1, const std::string &from2,
                                   const std::string &to, const double rate,
                                   const double eductDistance, const double weight1, const double weight2) {
    ensureKernelSelected();
    auto reaction = pimpl->kernel->createFusionReaction(name, from1, from2, to, rate, eductDistance, weight1, weight2);
    return pimpl->kernel->getKernelContext().registerReaction(std::move(reaction));
}

const short
Simulation::registerDecayReaction(const std::string &name, const std::string &particleType, const double rate) {
    ensureKernelSelected();
    auto reaction = pimpl->kernel->createDecayReaction(name, particleType, rate);
    return pimpl->kernel->getKernelContext().registerReaction(std::move(reaction));
}

std::vector<readdy::model::Vec3> Simulation::getParticlePositions(std::string type) {
    unsigned int typeId = pimpl->kernel->getKernelContext().getParticleTypeID(type);
    const auto particles = pimpl->kernel->getKernelStateModel().getParticles();
    std::vector<readdy::model::Vec3> positions;
    for (auto &&p : particles) {
        if (p.getType() == typeId) {
            positions.push_back(p.getPos());
        }
    }
    return positions;
}

double Simulation::getRecommendedTimeStep(unsigned int N) const {
    return readdy::model::util::getRecommendedTimeStep(N, pimpl->kernel->getKernelContext());
}

readdy::model::Kernel *const Simulation::getSelectedKernel() const {
    return pimpl->kernel.get();
}

void
Simulation::recordTrajectory(const std::string &fileName, const unsigned int stride, const unsigned int flushStride) {
    ensureKernelSelected();
    auto uuid = pimpl->counter++;
    pimpl->trajectoryFileId = uuid;
    pimpl->trajectoryFile.reset(new io::File(fileName, io::File::Action::CREATE));
    std::unique_ptr<model::observables::Trajectory> trajectory = std::make_unique<model::observables::Trajectory>(
            pimpl->kernel.get(), stride, flushStride, *pimpl->trajectoryFile
    );
    auto &&connection = pimpl->kernel->connectObservable(trajectory.get());
    pimpl->observables.emplace(uuid, std::move(trajectory));
    pimpl->observableConnections.emplace(uuid, std::move(connection));
}

void Simulation::closeTrajectoryFile() {
    deregisterObservable(pimpl->trajectoryFileId);
    pimpl->trajectoryFile.reset();
}

NoKernelSelectedException::NoKernelSelectedException(const std::string &__arg) : runtime_error(__arg) {};

}