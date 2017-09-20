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
#include <readdy/common/signals.h>
#include <readdy/model/Plugin.h>
#include <readdy/model/actions/Action.h>
#include <readdy/model/KernelStateModel.h>
#include <readdy/model/KernelContext.h>
#include <readdy/model/observables/ObservableFactory.h>
#include <readdy/model/_internal/ObservableWrapper.h>
#include <readdy/model/actions/ActionFactory.h>
#include <readdy/model/topologies/TopologyActionFactory.h>
#include <readdy/model/_internal/Util.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

NAMESPACE_BEGIN(detail)
template<typename T, typename... Args>
struct get_reaction_dispatcher;
NAMESPACE_END(detail)

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
    explicit Kernel(const std::string &name);

    /**
     * The kernel destructor.
     */
    ~Kernel() override;

    Kernel(const Kernel &rhs) = delete;

    Kernel &operator=(const Kernel &rhs) = delete;

    Kernel(Kernel &&rhs) noexcept;

    Kernel &operator=(Kernel &&rhs) noexcept;

    /**
     * This method returns the name of the kernel.
     *
     * @return The name
     */
    const std::string &getName() const override;

    template<typename ActionType, typename... Args>
    std::unique_ptr<ActionType> createAction(Args &&... args) const {
        return getActionFactory().createAction<ActionType>(std::forward<Args>(args)...);
    }

    template<typename T, typename... Args>
    std::unique_ptr<T> createObservable(unsigned int stride, Args &&... args) {
        return getObservableFactory().create<T>(stride, std::forward<Args>(args)...);
    }

    template<typename T, typename Obs1, typename Obs2>
    std::unique_ptr<T> createObservable(Obs1 *obs1, Obs2 *obs2, unsigned int stride = 1) {
        return getObservableFactory().create<T>(obs1, obs2, stride);
    }

    /**
     * Connects an observable to the kernel signal.
     *
     * @return A connection object that, once deleted, releases the connection of the observable.
     */
    virtual readdy::signals::scoped_connection connectObservable(observables::ObservableBase *observable);

    /**
     * Evaluates all observables.
     */
    virtual void evaluateObservables(time_step_type t);

    /**
     * Registers an observable to the kernel signal.
     */
    virtual std::tuple<std::unique_ptr<observables::ObservableWrapper>, readdy::signals::scoped_connection>
    registerObservable(const observables::observable_type &observable, unsigned int stride);


    /**
     * Returns a vector containing all available program names for this specific kernel instance.
     *
     * @see createProgram(name)
     * @return The program names.
     */
    virtual std::vector<std::string> getAvailableActions() const {
        return getActionFactory().getAvailableActions();
    }

    /**
     * Adds a particle of the type "type" at position "pos".
     */
    readdy::model::Particle::id_type addParticle(const std::string &type, const Vec3 &pos);

    TopologyParticle createTopologyParticle(const std::string &type, const Vec3 &pos) const;

    template<typename T, typename Obs1, typename Obs2>
    inline std::unique_ptr<T> createCombinedObservable(Obs1 *obs1, Obs2 *obs2, unsigned int stride = 1) const {
        return getObservableFactory().create(obs1, obs2, stride);
    }

    template<typename R, typename... Args>
    inline std::unique_ptr<R> createObservable(unsigned int stride, Args &&... args) const {
        return getObservableFactory().create(stride, std::forward<Args>(args)...);
    }

    bool supportsTopologies() const;

    /*
     * 
     * Accessors
     * 
     */

    const readdy::model::KernelStateModel &getKernelStateModel() const;

    readdy::model::KernelStateModel &getKernelStateModel();

    const readdy::model::KernelContext &getKernelContext() const;

    readdy::model::KernelContext &getKernelContext();

    const readdy::model::actions::ActionFactory &getActionFactory() const;

    readdy::model::actions::ActionFactory &getActionFactory();

    const readdy::model::observables::ObservableFactory &getObservableFactory() const;

    readdy::model::observables::ObservableFactory &getObservableFactory();

    const readdy::model::top::TopologyActionFactory *const getTopologyActionFactory() const;

    readdy::model::top::TopologyActionFactory *const getTopologyActionFactory();

    virtual void initialize();

    virtual void finalize();

    bool singlePrecision() const noexcept {
        return readdy::single_precision;
    }

    bool doublePrecision() const noexcept {
        return readdy::double_precision;
    }

protected:

    particle_type_type getTypeId(const std::string& name) const;

    particle_type_type getTypeIdRequireNormalFlavor(const std::string &) const;

    virtual readdy::model::KernelStateModel &getKernelStateModelInternal() const = 0;

    virtual readdy::model::KernelContext &getKernelContextInternal() const = 0;

    virtual readdy::model::actions::ActionFactory &getActionFactoryInternal() const = 0;

    virtual readdy::model::observables::ObservableFactory &getObservableFactoryInternal() const;

    virtual readdy::model::top::TopologyActionFactory *getTopologyActionFactoryInternal() const = 0;


    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
