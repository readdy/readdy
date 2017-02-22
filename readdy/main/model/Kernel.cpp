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


/**
 * @file Kernel.cpp
 * @brief Core library Implementation of the kernel.
 * @author clonker
 * @date 02.05.16
 */

#include <readdy/common/make_unique.h>
#include <readdy/model/Kernel.h>
#include <atomic>

namespace readdy {
namespace model {

struct Kernel::Impl {
    /**
     * The name of the kernel.
     */
    std::string name;
    /**
     * todo
     */
    std::unique_ptr<observables::signal_type> signal;
    /**
     * todo
     */
    std::unique_ptr<observables::ObservableFactory> observableFactory;

};

const std::string &Kernel::getName() const {
    return pimpl->name;
}

Kernel::Kernel(const std::string &name) : pimpl(std::make_unique<Kernel::Impl>()) {
    pimpl->name = name;
    pimpl->observableFactory = std::make_unique<observables::ObservableFactory>(this);
    pimpl->signal = std::make_unique<observables::signal_type>();
}

Kernel::~Kernel() {
}

readdy::signals::scoped_connection Kernel::connectObservable(observables::ObservableBase *const observable) {
    return pimpl->signal->connect_scoped([observable](const observables::time_step_type t) {
        observable->callback(t);
    });
}

std::tuple<std::unique_ptr<readdy::model::observables::ObservableWrapper>, readdy::signals::scoped_connection>
Kernel::registerObservable(const observables::observable_type &observable, unsigned int stride) {
    auto &&wrap = std::make_unique<observables::ObservableWrapper>(this, observable, stride);
    auto &&connection = connectObservable(wrap.get());
    return std::make_tuple(std::move(wrap), std::move(connection));
}

void Kernel::evaluateObservables(observables::time_step_type t) {
    (*pimpl->signal)(t);
}

std::vector<std::string> Kernel::getAvailablePotentials() const {
    return std::vector<std::string>();
}

void Kernel::addParticle(const std::string &type, const Vec3 &pos) {
    getKernelStateModel().addParticle({pos[0], pos[1], pos[2], getKernelContext().getParticleTypeID(type)});
}

unsigned int Kernel::getTypeIdRequireNormalFlavor(const std::string &name) const {
    auto findIt = getKernelContext().getTypeMapping().find(name);
    if (findIt != getKernelContext().getTypeMapping().end()) {
        const auto &info = getKernelContext().getParticleTypeInfo(findIt->second);
        if (info.flavor == readdy::model::Particle::FLAVOR_NORMAL) {
            return findIt->second;
        } else {
            log::critical("particle type {} had no \"normal\" flavor", name);
            throw std::invalid_argument("particle type " + name + " had no \"normal\" flavor");
        }
    } else {
        log::critical("did not find type id for {}", name);
        throw std::invalid_argument("did not find type id for " + name);
    }
}

std::unique_ptr<reactions::Reaction<1>>
Kernel::createConversionReaction(const std::string &name, const std::string &from, const std::string &to,
                                 const double rate) const {
    const auto typeFrom = getTypeIdRequireNormalFlavor(from);
    const auto typeTo = getTypeIdRequireNormalFlavor(to);
    return getReactionFactory().createReaction<reactions::Conversion>(name, typeFrom, typeTo, rate);
}

std::unique_ptr<reactions::Reaction<2>>
Kernel::createFusionReaction(const std::string &name, const std::string &from1, const std::string &from2,
                             const std::string &to, const double rate, const double eductDistance,
                             const double weight1, const double weight2) const {
    return getReactionFactory().createReaction<reactions::Fusion>(name, getTypeIdRequireNormalFlavor(from1),
                                                                  getTypeIdRequireNormalFlavor(from2),
                                                                  getTypeIdRequireNormalFlavor(to), rate, eductDistance,
                                                                  weight1, weight2);
}

std::unique_ptr<reactions::Reaction<2>>
Kernel::createEnzymaticReaction(const std::string &name, const std::string &catalyst, const std::string &from,
                                const std::string &to, const double rate, const double eductDistance) const {
    return getReactionFactory().createReaction<reactions::Enzymatic>(name, getTypeIdRequireNormalFlavor(catalyst),
                                                                     getTypeIdRequireNormalFlavor(from),
                                                                     getTypeIdRequireNormalFlavor(to), rate,
                                                                     eductDistance);
}

std::unique_ptr<reactions::Reaction<1>>
Kernel::createFissionReaction(const std::string &name, const std::string &from, const std::string &to1,
                              const std::string &to2, const double rate, const double productDistance,
                              const double weight1, const double weight2) const {
    return getReactionFactory().createReaction<reactions::Fission>(name, getTypeIdRequireNormalFlavor(from),
                                                                   getTypeIdRequireNormalFlavor(to1),
                                                                   getTypeIdRequireNormalFlavor(to2), rate,
                                                                   productDistance, weight1,
                                                                   weight2);
}

std::unique_ptr<reactions::Reaction<1>>
Kernel::createDecayReaction(const std::string &name, const std::string &type, const double rate) const {
    return getReactionFactory().createReaction<reactions::Decay>(name, getTypeIdRequireNormalFlavor(type), rate);
}

observables::ObservableFactory &Kernel::getObservableFactoryInternal() const {
    return *pimpl->observableFactory;
}

const readdy::model::KernelContext &Kernel::getKernelContext() const {
    return getKernelContextInternal();
}

readdy::model::KernelContext &Kernel::getKernelContext() {
    return getKernelContextInternal();
}

const readdy::model::KernelStateModel &Kernel::getKernelStateModel() const {
    return getKernelStateModelInternal();
}

readdy::model::KernelStateModel &Kernel::getKernelStateModel() {
    return getKernelStateModelInternal();
}

const readdy::model::actions::ActionFactory &Kernel::getActionFactory() const {
    return getActionFactoryInternal();
}

readdy::model::actions::ActionFactory &Kernel::getActionFactory() {
    return getActionFactoryInternal();
}

const readdy::model::potentials::PotentialFactory &Kernel::getPotentialFactory() const {
    return getPotentialFactoryInternal();
}

readdy::model::potentials::PotentialFactory &Kernel::getPotentialFactory() {
    return getPotentialFactoryInternal();
}

const readdy::model::reactions::ReactionFactory &Kernel::getReactionFactory() const {
    return getReactionFactoryInternal();
}

readdy::model::reactions::ReactionFactory &Kernel::getReactionFactory() {
    return getReactionFactoryInternal();
}

const readdy::model::compartments::CompartmentFactory &Kernel::getCompartmentFactory() const {
    return getCompartmentFactoryInternal();
}

readdy::model::compartments::CompartmentFactory &Kernel::getCompartmentFactory() {
    return getCompartmentFactoryInternal();
}

const readdy::model::observables::ObservableFactory &Kernel::getObservableFactory() const {
    return getObservableFactoryInternal();
}

readdy::model::observables::ObservableFactory &Kernel::getObservableFactory() {
    return getObservableFactoryInternal();
}

const readdy::model::top::TopologyActionFactory *const Kernel::getTopologyActionFactory() const {
    return getTopologyActionFactoryInternal();
}

readdy::model::top::TopologyActionFactory *const Kernel::getTopologyActionFactory() {
    return getTopologyActionFactoryInternal();
}

TopologyParticle Kernel::createTopologyParticle(const std::string &type, const Vec3 &pos) const {
    return TopologyParticle(pos, getKernelContext().getParticleTypeID(type));
}

bool Kernel::supportsTopologies() const {
    return getTopologyActionFactory() ? true : false;
}

Kernel &Kernel::operator=(Kernel &&rhs) = default;

Kernel::Kernel(Kernel &&rhs) = default;
}
}
