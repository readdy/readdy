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
 * A Kernel is used to execute Programs, i.e., instances of readdy::plugin::Program.
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

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

/**
 * Base class of kernels.
 * A Kernel is used to execute Actions, i.e., instances of readdy::model::actions::Action.
 * The kernels can be built in or provided by shared libs in directories, which are loaded by the KernelProvider.
 * Each Kernel has a #name by which it can be accessed in the KernelProvider.
 */
class Kernel : public Plugin {
public:

    /**
     * Constructs a kernel with a given name.
     */
    explicit Kernel(std::string name) :_name(std::move(name)), _signal() {};

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
    };

    /**
     * Connects an observable to the kernel signal.
     *
     * @return A connection object that, once deleted, releases the connection of the observable.
     */
    virtual readdy::signals::scoped_connection connectObservable(observables::ObservableBase *observable) {
        observable->initialize(this);
        return _signal.connect_scoped([observable](const time_step_type t) {
            observable->callback(t);
        });
    };

    /**
     * Evaluates all observables.
     */
    virtual void evaluateObservables(time_step_type t) {
        _signal(t);
    };

    /**
     * Returns a vector containing all available program names for this specific kernel instance.
     *
     * @see createProgram(name)
     * @return The program names.
     */
    virtual std::vector<std::string> getAvailableActions() const {
        return actions().getAvailableActions();
    }

    /**
     * Adds a particle of the type "type" at position "pos".
     */
    readdy::model::Particle::id_type addParticle(const std::string &type, const Vec3 &pos) {
        readdy::model::Particle particle {pos[0], pos[1], pos[2], context().particleTypes().idOf(type)};
        stateModel().addParticle(particle);
        return particle.id();
    };

    TopologyParticle createTopologyParticle(const std::string &type, const Vec3 &pos) const {
        const auto& info = context().particleTypes().infoOf(type);
        if(info.flavor != particleflavor::TOPOLOGY) {
            throw std::invalid_argument("You can only create topology particles of a type that is topology flavored.");
        }
        return TopologyParticle(pos, info.typeId);
    };

    bool supportsTopologies() const {
        return getTopologyActionFactory() != nullptr;
    };

    /*
     * 
     * Accessors
     * 
     */

    virtual const readdy::model::StateModel &stateModel() const = 0;

    virtual readdy::model::StateModel &stateModel() = 0;

    const readdy::model::Context &context() const {
        return _context;
    };

    readdy::model::Context &context() {
        return _context;
    };

    virtual const readdy::model::actions::ActionFactory &actions() const = 0;

    virtual readdy::model::actions::ActionFactory &actions() = 0;

    virtual const readdy::model::observables::ObservableFactory &observe() const = 0;

    virtual readdy::model::observables::ObservableFactory &observe() = 0;

    virtual const readdy::model::top::TopologyActionFactory *const getTopologyActionFactory() const = 0;

    virtual readdy::model::top::TopologyActionFactory *const getTopologyActionFactory() = 0;

    virtual void initialize() {
        log::debug(context().describe());
    };

    bool singlePrecision() const noexcept {
        return readdy::single_precision;
    }

    bool doublePrecision() const noexcept {
        return readdy::double_precision;
    }

protected:

    model::Context _context;
    std::string _name;
    observables::signal_type _signal;
};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
