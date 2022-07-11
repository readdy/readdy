/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * This file contains the base class definitions for Kernel and KernelProvider.
 * A Kernel is used to execute Actions, i.e., instances of readdy::model::Action.
 * The derived kernels can be built in or provided by shared libs in directories, which are loaded by the KernelProvider.
 * Each Kernel has a readdy::plugin::Kernel::name by which it can be accessed in the KernelProvider.
 *
 * @file Kernel.h
 * @brief File containing the definitions of Kernel and KernelProvider.
 * @author clonker
 * @date 07.03.16
 */

#pragma once

#include <map>
#include <iostream>
#include <utility>
#include <readdy/common/signals.h>
#include <readdy/model/Plugin.h>
#include <readdy/model/actions/Action.h>
#include <readdy/model/StateModel.h>
#include <readdy/model/Context.h>
#include <readdy/model/observables/ObservableFactory.h>
#include <readdy/model/actions/ActionFactory.h>
#include <readdy/model/topologies/TopologyActionFactory.h>
#include <readdy/model/_internal/Util.h>

namespace readdy::model {

/**
 * Base class of kernels.
 * A Kernel is used to execute Actions, i.e., instances of readdy::model::actions::Action.
 * The kernels can be built in or provided by shared libs in directories, which are loaded by the KernelProvider.
 * Each Kernel has a #name by which it can be accessed in the KernelProvider.
 */
class Kernel : public Plugin {
    using ObservableContainer = std::vector<std::unique_ptr<readdy::model::observables::ObservableBase>>;
    using ConnectionContainer = std::vector<readdy::signals::scoped_connection>;
public:

    /**
     * Constructs a kernel with a given name.
     */
    explicit Kernel(std::string name) : Kernel(std::move(name), {}) {}

    Kernel(std::string name, readdy::model::Context ctx) : _name(std::move(name)), _signal(), _context(std::move(ctx)) {}

    /**
     * The kernel destructor.
     */
    ~Kernel() override = default;

    Kernel(const Kernel &rhs) = delete;

    Kernel &operator=(const Kernel &rhs) = delete;

    /**
     * This method returns the name of the kernel.
     *
     * @return The name
     */
    const std::string &name() const override {
        return _name;
    }

    /**
     * Connect the observable, store the connection and the observable. As long as the connection stays alive,
     * the observable will be evaluated upon firing the signal. Return a handle to the observable, e.g. used for
     * saving results to file.
     * @param observable
     * @return an observable handle that allows for post-hoc modification of the observable
     */
    ObservableHandle registerObservable(std::unique_ptr<readdy::model::observables::ObservableBase> observable) {
        auto connection = connectObservable(observable.get());
        registeredObservables().push_back(std::move(observable));
        observableConnections().push_back(std::move(connection));
        return ObservableHandle{registeredObservables().back().get()};
    }

    /**
     * Connects an observable to the kernel signal.
     *
     * @return A connection object that, once deleted, releases the connection of the observable.
     */
    virtual readdy::signals::scoped_connection connectObservable(observables::ObservableBase *observable) {
        observable->initialize(this);
        return _signal.connect_scoped([observable](const TimeStep t) {
            observable->call(t);
        });
    }

    /**
     * Evaluates all observables.
     */
    virtual void evaluateObservables(TimeStep t) {
        _signal(t);
    }

    /**
     * Returns a vector containing all available action names for this specific kernel instance.
     *
     * @return The action names.
     */
    virtual std::vector<std::string> getAvailableActions() const {
        return actions().getAvailableActions();
    }

    /**
     * Adds a particle of the type "type" at position "pos".
     */
    ParticleId addParticle(const std::string &type, const Vec3 &pos) {
        readdy::model::Particle particle{pos[0], pos[1], pos[2], context().particleTypes().idOf(type)};
        stateModel().addParticle(particle);
        return particle.id();
    }

    Particle createTopologyParticle(const std::string &type, const Vec3 &pos) const {
        const auto &info = context().particleTypes().infoOf(type);
        if (info.flavor != particleflavor::TOPOLOGY) {
            throw std::invalid_argument(fmt::format(
                    "You can only create topology particles of a type that is topology flavored. "
                    "Type was {}, flavor {}.", type, readdy::model::particleflavor::particleFlavorToString(info.flavor)
            ));
        }
        return Particle(pos, info.typeId);
    }

    bool supportsTopologies() const {
        return getTopologyActionFactory() != nullptr;
    }
    
    virtual bool supportsGillespie() const {
        return true;
    }

    /*
     * 
     * Accessors
     * 
     */

    virtual const readdy::model::StateModel &stateModel() const = 0;

    virtual readdy::model::StateModel &stateModel() = 0;

    const readdy::model::Context &context() const {
        return _context;
    }

    readdy::model::Context &context() {
        return _context;
    }

    virtual const readdy::model::actions::ActionFactory &actions() const = 0;

    virtual readdy::model::actions::ActionFactory &actions() = 0;

    virtual const readdy::model::observables::ObservableFactory &observe() const = 0;

    virtual readdy::model::observables::ObservableFactory &observe() = 0;

    virtual const readdy::model::top::TopologyActionFactory *const getTopologyActionFactory() const = 0;

    virtual readdy::model::top::TopologyActionFactory *const getTopologyActionFactory() = 0;

    virtual void initialize() {
        log::debug("{}", context().describe());
    }

    bool singlePrecision() const noexcept {
        return readdy::single_precision;
    }

    bool doublePrecision() const noexcept {
        return readdy::double_precision;
    }

    const std::vector<std::unique_ptr<readdy::model::observables::ObservableBase>> &registeredObservables() const {
        return _observables;
    }

    std::vector<std::unique_ptr<readdy::model::observables::ObservableBase>> &registeredObservables() {
        return _observables;
    }

    const std::vector<readdy::signals::scoped_connection> &observableConnections() const {
        return _observableConnections;
    }

    std::vector<readdy::signals::scoped_connection> &observableConnections() {
        return _observableConnections;
    }

protected:
    model::Context _context;
    std::string _name;
    observables::signal_type _signal;
    ObservableContainer _observables{};
    ConnectionContainer _observableConnections{};
};

}
