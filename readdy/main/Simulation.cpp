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
#include <readdy/model/Utils.h>
#include <readdy/model/observables/io/Trajectory.h>

namespace readdy {
scalar Simulation::getKBT() const {
    ensureKernelSelected();
    return pimpl->kernel->getKernelContext().getKBT();

}

void Simulation::setKBT(scalar kBT) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().setKBT(kBT);

}

void Simulation::setBoxSize(const readdy::model::Vec3 &boxSize) {
    setBoxSize(boxSize[0], boxSize[1], boxSize[2]);
}

void Simulation::setPeriodicBoundary(const std::array<bool, 3> &periodic) {
    ensureKernelSelected();
    pimpl->kernel->getKernelContext().setPeriodicBoundary(periodic[0], periodic[1], periodic[2]);

}

Simulation::Simulation() : pimpl(std::make_unique<Simulation::Impl>()) {}

readdy::model::Vec3 Simulation::getBoxSize() const {
    ensureKernelSelected();
    return readdy::model::Vec3(pimpl->kernel->getKernelContext().getBoxSize());

}

const std::array<bool, 3>& Simulation::getPeriodicBoundary() const {
    ensureKernelSelected();
    return pimpl->kernel->getKernelContext().getPeriodicBoundary();
}

void Simulation::run(const time_step_type steps, const scalar timeStep) {
    ensureKernelSelected();
    {
        log::debug("available actions: ");
        for (auto &&p : pimpl->kernel->getAvailableActions()) {
            log::debug("\t {}", p);
        }
    }
    runScheme().configure(timeStep)->run(steps);
}

void Simulation::setKernel(const std::string &kernel) {
    if (isKernelSelected()) {
        log::warn(R"(replacing kernel "{}" with "{}")", pimpl->kernel->getName(), kernel);
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

void Simulation::addParticle(const std::string &type, scalar x, scalar y, scalar z) {
    ensureKernelSelected();
    const auto &&s = getBoxSize();
    if (fabs(x) <= .5 * s[0] && fabs(y) <= .5 * s[1] && fabs(z) <= .5 * s[2]) {
        readdy::model::Particle p{x, y, z, pimpl->kernel->getKernelContext().particle_types().id_of(type)};
        pimpl->kernel->getKernelStateModel().addParticle(p);
    } else {
        log::error("particle position was not in bounds of the simulation box!");
    }

}

Simulation::particle::type_type
Simulation::registerParticleType(const std::string &name, const scalar diffusionCoefficient, const scalar radius,
                                 readdy::model::particle_flavor flavor) {
    ensureKernelSelected();
    auto &context = pimpl->kernel->getKernelContext();
    context.particle_types().add(name, diffusionCoefficient, radius, flavor);
    return context.particle_types().id_of(name);
}

const std::vector<readdy::model::Vec3> Simulation::getAllParticlePositions() const {
    ensureKernelSelected();
    return pimpl->kernel->getKernelStateModel().getParticlePositions();
}

void Simulation::deregisterPotential(const short uuid) {
    pimpl->kernel->getKernelContext().potentials().remove(uuid);
}

const short
Simulation::registerHarmonicRepulsionPotential(const std::string &particleTypeA, const std::string &particleTypeB,
                                               scalar forceConstant) {
    using potential_t = readdy::model::potentials::HarmonicRepulsion;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleTypeA, particleTypeB, forceConstant);
}

const short
Simulation::registerWeakInteractionPiecewiseHarmonicPotential(const std::string &particleTypeA,
                                                              const std::string &particleTypeB, scalar forceConstant,
                                                              scalar desiredParticleDistance, scalar depth,
                                                              scalar noInteractionDistance) {
    using potential_t = readdy::model::potentials::WeakInteractionPiecewiseHarmonic;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleTypeA, particleTypeB, forceConstant,
                                                         desiredParticleDistance, depth, noInteractionDistance);
}

const short
Simulation::registerLennardJonesPotential(const std::string &type1, const std::string &type2, unsigned int m,
                                          unsigned int n, scalar cutoff, bool shift, scalar epsilon, scalar sigma) {
    using potential_t = readdy::model::potentials::LennardJones;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(type1, type2, m, n, cutoff, shift, epsilon, sigma);
}


const short
Simulation::registerScreenedElectrostaticsPotential(const std::string &particleType1, const std::string &particleType2,
                                                    scalar electrostaticStrength,
                                                    scalar inverseScreeningDepth, scalar repulsionStrength,
                                                    scalar repulsionDistance,
                                                    unsigned int exponent, scalar cutoff) {
    using potential_t = readdy::model::potentials::ScreenedElectrostatics;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleType1, particleType2, electrostaticStrength,
                                                         inverseScreeningDepth,
                                                         repulsionStrength, repulsionDistance, exponent, cutoff);
}

const short
Simulation::registerBoxPotential(const std::string &particleType, scalar forceConstant,
                                 const readdy::model::Vec3 &origin, const readdy::model::Vec3 &extent,
                                 bool considerParticleRadius) {
    using potential_t = readdy::model::potentials::Cube;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleType, forceConstant, origin, extent,
                                                         considerParticleRadius);
}

const short
Simulation::registerSphereInPotential(const std::string &particleType, scalar forceConstant, const readdy::model::Vec3 &origin,
                                      scalar radius) {
    using potential_t = readdy::model::potentials::SphereIn;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleType, forceConstant, origin, radius);
}

const short
Simulation::registerSphereOutPotential(const std::string &particleType, scalar forceConstant,
                                       const readdy::model::Vec3 &origin, scalar radius) {
    using potential_t = readdy::model::potentials::SphereOut;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleType, forceConstant, origin, radius);
}

const short
Simulation::registerSphericalBarrier(const std::string &particleType, const readdy::model::Vec3 &origin, scalar radius, scalar height, scalar width) {
    using potential_t = readdy::model::potentials::SphericalBarrier;
    ensureKernelSelected();
    return pimpl->kernel->registerPotential<potential_t>(particleType, origin, radius, height, width);
}

void Simulation::ensureKernelSelected() const {
    if (!isKernelSelected()) {
        throw NoKernelSelectedException("No kernel was selected!");
    }
}

void Simulation::setBoxSize(scalar dx, scalar dy, scalar dz) {
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

Simulation &Simulation::operator=(Simulation &&rhs) noexcept = default;

Simulation::Simulation(Simulation &&rhs) noexcept = default;

Simulation::~Simulation() {
    log::trace("destroying simulation");
}

const short
Simulation::registerConversionReaction(const std::string &name, const std::string &from, const std::string &to,
                                       const scalar rate) {
    ensureKernelSelected();
    auto reaction = pimpl->kernel->createConversionReaction(name, from, to, rate);
    return pimpl->kernel->getKernelContext().reactions().add(std::move(reaction));
}

const short
Simulation::registerEnzymaticReaction(const std::string &name, const std::string &catalyst, const std::string &from,
                                      const std::string &to, const scalar rate,
                                      const scalar eductDistance) {
    ensureKernelSelected();
    auto reaction = pimpl->kernel->createEnzymaticReaction(name, catalyst, from, to, rate, eductDistance);
    return pimpl->kernel->getKernelContext().reactions().add(std::move(reaction));
}

const short
Simulation::registerFissionReaction(const std::string &name, const std::string &from, const std::string &to1,
                                    const std::string &to2,
                                    const scalar rate, const scalar productDistance, const scalar weight1,
                                    const scalar weight2) {
    ensureKernelSelected();
    auto reaction = pimpl->kernel->createFissionReaction(name, from, to1, to2, rate, productDistance, weight1, weight2);
    return pimpl->kernel->getKernelContext().reactions().add(std::move(reaction));
}

const short
Simulation::registerFusionReaction(const std::string &name, const std::string &from1, const std::string &from2,
                                   const std::string &to, const scalar rate,
                                   const scalar eductDistance, const scalar weight1, const scalar weight2) {
    ensureKernelSelected();
    auto reaction = pimpl->kernel->createFusionReaction(name, from1, from2, to, rate, eductDistance, weight1, weight2);
    return pimpl->kernel->getKernelContext().reactions().add(std::move(reaction));
}

const short
Simulation::registerDecayReaction(const std::string &name, const std::string &particleType, const scalar rate) {
    ensureKernelSelected();
    auto reaction = pimpl->kernel->createDecayReaction(name, particleType, rate);
    return pimpl->kernel->getKernelContext().reactions().add(std::move(reaction));
}

std::vector<readdy::model::Vec3> Simulation::getParticlePositions(std::string type) {
    unsigned int typeId = pimpl->kernel->getKernelContext().particle_types().id_of(type);
    const auto particles = pimpl->kernel->getKernelStateModel().getParticles();
    std::vector<readdy::model::Vec3> positions;
    for (auto &&p : particles) {
        if (p.getType() == typeId) {
            positions.push_back(p.getPos());
        }
    }
    return positions;
}

scalar Simulation::getRecommendedTimeStep(unsigned int N) const {
    return readdy::model::util::getRecommendedTimeStep(N, pimpl->kernel->getKernelContext());
}

readdy::model::Kernel *const Simulation::getSelectedKernel() const {
    ensureKernelSelected();
    return pimpl->kernel.get();
}

const short Simulation::registerCompartmentSphere(const std::unordered_map<std::string, std::string> &conversionsMap,
                                                  const std::string &name,
                                                  const model::Vec3 &origin, const scalar radius,
                                                  const bool largerOrLess) {
    ensureKernelSelected();
    return getSelectedKernel()->registerCompartment<model::compartments::Sphere>(conversionsMap, name, origin, radius,
                                                                                 largerOrLess);
}

const short Simulation::registerCompartmentPlane(const std::unordered_map<std::string, std::string> &conversionsMap,
                                                 const std::string &name,
                                                 const model::Vec3 &normalCoefficients, const scalar distanceFromPlane,
                                                 const bool largerOrLess) {
    ensureKernelSelected();
    return getSelectedKernel()->registerCompartment<model::compartments::Plane>(conversionsMap, name,
                                                                                normalCoefficients, distanceFromPlane,
                                                                                largerOrLess);
}

readdy::model::TopologyParticle
Simulation::createTopologyParticle(const std::string &type, const readdy::model::Vec3 &pos) const {
    ensureKernelSelected();
    return getSelectedKernel()->createTopologyParticle(type, pos);
}

bool Simulation::kernelSupportsTopologies() const {
    ensureKernelSelected();
    return getSelectedKernel()->supportsTopologies();
}

readdy::model::top::GraphTopology *
Simulation::addTopology(const std::string& type, const std::vector<readdy::model::TopologyParticle> &particles) {
    ensureKernelSelected();
    if (getSelectedKernel()->supportsTopologies()) {
        auto typeId = getSelectedKernel()->getKernelContext().topology_registry().id_of(type);
        return getSelectedKernel()->getKernelStateModel().addTopology(typeId, particles);
    }
    throw std::logic_error("the selected kernel does not support topologies!");
}

void Simulation::registerPotentialOrder1(readdy::model::potentials::PotentialOrder1 *ptr) {
    ensureKernelSelected();
    getSelectedKernel()->getKernelContext().potentials().add_external(ptr);
}

void Simulation::registerPotentialOrder2(readdy::model::potentials::PotentialOrder2 *ptr) {
    ensureKernelSelected();
    getSelectedKernel()->getKernelContext().potentials().add_external(ptr);
}

void
Simulation::configureTopologyBondPotential(const std::string &type1, const std::string &type2, scalar forceConstant,
                                           scalar length, api::BondType type) {
    ensureKernelSelected();
    getSelectedKernel()->getKernelContext().topology_registry().configure_bond_potential(type1, type2,
                                                                                         {forceConstant, length, type});
}

void Simulation::configureTopologyAnglePotential(const std::string &type1, const std::string &type2,
                                                 const std::string &type3, scalar forceConstant,
                                                 scalar equilibriumAngle, api::AngleType type) {
    ensureKernelSelected();
    getSelectedKernel()->getKernelContext().topology_registry().configure_angle_potential(type1, type2, type3,
                                                                                          {forceConstant,
                                                                                           equilibriumAngle, type});
}

void Simulation::configureTopologyTorsionPotential(const std::string &type1, const std::string &type2,
                                                   const std::string &type3, const std::string &type4,
                                                   scalar forceConstant, unsigned int multiplicity, scalar phi_0,
                                                   api::TorsionType type) {
    ensureKernelSelected();
    getSelectedKernel()->getKernelContext().topology_registry().configure_torsion_potential(
            type1, type2, type3, type4, {forceConstant, static_cast<scalar>(multiplicity), phi_0, type}
    );
}

void Simulation::setExpectedMaxNParticles(const std::size_t n) {
    ensureKernelSelected();
    getSelectedKernel()->expected_n_particles(n);
}

std::vector<readdy::model::top::GraphTopology *> Simulation::currentTopologies() {
    ensureKernelSelected();
    return getSelectedKernel()->getKernelStateModel().getTopologies();
}

std::vector<model::Particle> Simulation::getParticlesForTopology(const model::top::GraphTopology &topology) const {
    ensureKernelSelected();
    return getSelectedKernel()->getKernelStateModel().getParticlesForTopology(topology);
}

bool Simulation::singlePrecision() const {
    ensureKernelSelected();
    return pimpl->kernel->singlePrecision();
}

bool Simulation::doublePrecision() const {
    ensureKernelSelected();
    return pimpl->kernel->doublePrecision();
}

void Simulation::registerSpatialTopologyReaction(const std::string &descriptor, scalar rate, scalar radius) {
    ensureKernelSelected();
    getSelectedKernel()->getKernelContext().topology_registry().add_spatial_reaction(descriptor, rate, radius);
}

readdy::plugin::KernelProvider::raw_kernel_ptr Simulation::setKernel(plugin::KernelProvider::kernel_ptr &&kernel) {
    if (isKernelSelected()) {
        log::warn(R"(replacing kernel "{}" with "{}")", pimpl->kernel->getName(), kernel->getName());
    }
    pimpl->kernel = std::move(kernel);
    return pimpl->kernel.get();
}

readdy::topology_type_type Simulation::registerTopologyType(const std::string &name,
                                      const std::vector<model::top::reactions::StructuralTopologyReaction> &reactions) {
    ensureKernelSelected();
    return getSelectedKernel()->getKernelContext().topology_registry().add_type(name, reactions);
}

void Simulation::registerStructuralTopologyReaction(const std::string &topologyType,
                                                    const model::top::reactions::StructuralTopologyReaction &reaction) {
    ensureKernelSelected();
    getSelectedKernel()->getKernelContext().topology_registry().add_structural_reaction(topologyType, reaction);
}

NoKernelSelectedException::NoKernelSelectedException(const std::string &__arg) : runtime_error(__arg) {}

}