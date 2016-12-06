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
 * << detailed description >>
 *
 * @file Kernel.cpp.cpp
 * @brief << brief description >>
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

readdy::model::observables::ObservableFactory &Kernel::getObservableFactory() const {
    return *pimpl->observableFactory;
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

std::unique_ptr<readdy::model::potentials::Potential> Kernel::createPotential(std::string &name) const {
    return getPotentialFactory().createPotential(name);
}

unsigned int Kernel::getTypeId(const std::string &name) const {
    return getKernelContext().getTypeMapping().find(name)->second;
}

std::unique_ptr<reactions::Reaction<1>>
Kernel::createConversionReaction(const std::string &name, const std::string &from, const std::string &to,
                                 const double rate) const {
    return getReactionFactory().createReaction<reactions::Conversion>(name, getTypeId(from), getTypeId(to), rate);
}

std::unique_ptr<reactions::Reaction<2>>
Kernel::createFusionReaction(const std::string &name, const std::string &from1, const std::string &from2,
                             const std::string &to, const double rate, const double eductDistance,
                             const double weight1, const double weight2) const {
    return getReactionFactory().createReaction<reactions::Fusion>(name, getTypeId(from1), getTypeId(from2),
                                                                  getTypeId(to), rate, eductDistance, weight1, weight2);
}

std::unique_ptr<reactions::Reaction<2>>
Kernel::createEnzymaticReaction(const std::string &name, const std::string &catalyst, const std::string &from,
                                const std::string &to, const double rate, const double eductDistance) const {
    return getReactionFactory().createReaction<reactions::Enzymatic>(name, getTypeId(catalyst), getTypeId(from),
                                                                     getTypeId(to), rate, eductDistance);
}

std::unique_ptr<reactions::Reaction<1>>
Kernel::createFissionReaction(const std::string &name, const std::string &from, const std::string &to1,
                              const std::string &to2, const double rate, const double productDistance,
                              const double weight1, const double weight2) const {
    return getReactionFactory().createReaction<reactions::Fission>(name, getTypeId(from), getTypeId(to1),
                                                                   getTypeId(to2), rate, productDistance, weight1,
                                                                   weight2);
}

std::unique_ptr<reactions::Reaction<1>>
Kernel::createDecayReaction(const std::string &name, const std::string &type, const double rate) const {
    return getReactionFactory().createReaction<reactions::Decay>(name, getTypeId(type), rate);
}


Kernel &Kernel::operator=(Kernel &&rhs) = default;

Kernel::Kernel(Kernel &&rhs) = default;
}
}



