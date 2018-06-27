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
 * Vital parts of the scheme api are defined in this header. First of all, there is the SimulationScheme. Its task is
 * to execute an integrator(action), a force(action), a reaction scheduler(action) and a neighborList (action)
 * in a certain way. One example would be the ReaDDyScheme, which will configure the context, then:
 *
 * - create neighbor list
 * - calculate forces
 * - evaluate observables (0)
 * - for i in range(steps):
 *      - evaluate integrator
 *      - update neighbor list
 *      - calculate forces
 *      - evaluate reactions
 *      - update the neighbor list
 *      - calculate forces
 *      - evaluate observables (i+1)
 *  - clear neighbor list
 *
 *  The configurator makes use of the builder pattern and either sets default values or not, depending on
 *  its constructor argument.
 *
 * @file SimulationScheme.h
 * @brief Header file containing the SchemeConfigurator<T> and a corresponding SimulationScheme superclass definition.
 * @author clonker
 * @author chrisfroe
 * @date 23.08.16
 */

#pragma once

#include <memory>
#include <type_traits>

#include <h5rd/h5rd.h>

#include <readdy/common/common.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/IOUtils.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(api)

/**
 * superclass for all simulation schemes
 */
class SimulationScheme {
public:
    /**
     * the type of function that is responsible for deciding whether the simulation should continue
     */
    using continue_fun = std::function<bool(time_step_type)>;

    /**
     * Creates a new simulation scheme.
     * @param kernel the kernel
     * @param performanceRoot performance stuff
     */
    explicit SimulationScheme(model::Kernel *const kernel, util::PerformanceNode &performanceRoot)
            : kernel(kernel), _performanceRoot(performanceRoot) {}

    /**
     * This method is the interface for running a simulation as long as a continue function yields true.
     * @param fun the continue function
     */
    virtual void run(const continue_fun &fun) = 0;

    /**
     * This method runs a simulation for a fixed number of time steps by providing an appropriate continue function.
     * @param steps the number of steps
     */
    virtual void run(time_step_type steps) {
        // show every 100 time steps
        if(!_updateCallback) {
            _updateCallback = [this, steps](time_step_type current) {
                log::info("Simulation progress: {} / {} steps", (current - start), steps);
            };
        }
        auto defaultContinueCriterion = [this, steps](const time_step_type current) {
            if (current != start && _progressOutputStride > 0 && (current - start) % _progressOutputStride == 0) {
                _updateCallback(current);
            }
            return current < start + steps;
        };
        run(defaultContinueCriterion);
    };

    /**
     * This function gives access the an updateCallback function that gets called every 100 time steps if one is
     * running the simulation for a fixed number of time steps.
     * @return the callback function
     */
    std::function<void(time_step_type)> &updateCallback() {
        return _updateCallback;
    }

    /**
     * show progress every N steps
     * @return N steps
     */
    std::size_t &progressOutputStride() {
        return _progressOutputStride;
    }



protected:
    friend class SchemeConfigurator;

    /**
     * reference to the kernel
     */
    model::Kernel *const kernel;
    /**
     * the integrator to use
     */
    std::shared_ptr<model::actions::TimeStepDependentAction> _integrator {nullptr};
    /**
     * how to evaluate forces
     */
    std::shared_ptr<model::actions::Action> _forces {nullptr};
    /**
     * how to evaluate particle-particle reactions
     */
    std::shared_ptr<model::actions::TimeStepDependentAction> _reactions {nullptr};
    /**
     * initialize the neighbor list
     */
    std::shared_ptr<model::actions::UpdateNeighborList> _initNeighborList {nullptr};
    /**
     * update the neighbor list
     */
    std::shared_ptr<model::actions::UpdateNeighborList> _neighborList {nullptr};
    /**
     * clear the neighbor list
     */
    std::shared_ptr<model::actions::UpdateNeighborList> _clearNeighborList {nullptr};
    /**
     * how to evaluate topology reactions (structural and spatial)
     */
    std::shared_ptr<model::actions::top::EvaluateTopologyReactions> _topologyReactions {nullptr};
    /**
     * reference to a hdf5-group holding some configurationn values
     */
    std::shared_ptr<h5rd::Group> configGroup = nullptr;
    /**
     * the update callback function, called every 100 time steps if running the simulation for a fixed number of steps
     */
    std::function<void(time_step_type)> _updateCallback;
    /**
     * wheter to evaluate observables at all, also affects writing trajectories
     */
    bool _evaluateObservables = true;
    /**
     * the starting point
     */
    time_step_type start = 0;
    /**
     * show progress every N steps
     */
    std::size_t _progressOutputStride = 100;
    /**
     * cref to the performance root node
     */
    const util::PerformanceNode &_performanceRoot;
};

/**
 * The ReaDDy scheme. By default:
 * - create neighbor list
 * - calculate forces
 * - evaluate observables (0)
 * - for i in range(steps):
 *      - evaluate integrator
 *      - update neighbor list
 *      - calculate forces
 *      - evaluate reactions
 *      - update the neighbor list
 *      - calculate forces
 *      - evaluate observables (i+1)
 *  - clear neighbor list
 */
class ReaDDyScheme : public SimulationScheme {
public:
    /**
     * Instantiates a new ReaDDy scheme
     * @param kernel the kernel
     * @param performanceRoot perf
     */
    explicit ReaDDyScheme(model::Kernel *const kernel, util::PerformanceNode &performanceRoot)
            : SimulationScheme(kernel, performanceRoot) {};

    using SimulationScheme::run;

    void initialize() {
        kernel->initialize();
        if(configGroup) {
            model::ioutils::writeSimulationSetup(*configGroup, kernel->context());
        }
    }

    void initializeNeighborList() {
        if (_initNeighborList) _initNeighborList->perform(_performanceRoot.subnode("initNeighborList"));
    }

    void updateNeighborList() {
        if (_neighborList) _neighborList->perform(_performanceRoot.subnode("neighborList"));
    }

    void clearNeighborList() {
        if (_clearNeighborList) _clearNeighborList->perform(_performanceRoot.subnode("clearNeighborList"));
    }

    void forces() {
        if (_forces) _forces->perform(_performanceRoot.subnode("forces"));
    }

    void evaluateObservables(time_step_type t) {
        if (_evaluateObservables) kernel->evaluateObservables(t);
    }

    void integrator() {
        if (_integrator) _integrator->perform(_performanceRoot.subnode("integrator"));
    }

    void reactions() {
        if (_reactions) _reactions->perform(_performanceRoot.subnode("reactionScheduler"));
    }

    void topologyReactions() {
        if (_topologyReactions) _topologyReactions->perform(_performanceRoot.subnode("evaluateTopologyReactions"));
    }

    /**
     * ReaDDy scheme implementation of the simulation loop
     * @param continueFun the continue function
     */
    void run(const continue_fun &continueFun) override {
        auto runTimer = _performanceRoot.timeit();
        initialize();
        initializeNeighborList();
        forces();
        time_step_type t = start;
        evaluateObservables(t);
        while (continueFun(t)) {
            integrator();
            updateNeighborList();
            reactions();
            topologyReactions();
            updateNeighborList();
            forces();
            evaluateObservables(t+1);
            ++t;
        }
        clearNeighborList();
        start = t;
        log::info("Simulation completed");
    }
};

namespace detail {
template<typename T>
struct identity { typedef T type; };
}

/**
 * The scheme configurator, enabling fluent-api style configuration of the simulation scheme.
 * @tparam SchemeType the scheme type
 */
class SchemeConfigurator {
public:

    /**
     * Instantiates a new scheme configurator object by providing a kernel, the performance root and whether to use
     * default values.
     * @param kernel the kernel
     * @param perfRoot perf
     * @param useDefaults whether to use defaults
     */
    explicit SchemeConfigurator(model::Kernel *const kernel, util::PerformanceNode &perfRoot, bool useDefaults = true)
            : scheme(std::make_unique<ReaDDyScheme>(kernel, perfRoot)), useDefaults(useDefaults) {}

    /**
     * Instruct the configurator that the simulation loop should use the provided integrator instance
     * @param integrator the integrator
     * @return reference to self
     */
    SchemeConfigurator &withIntegrator(std::shared_ptr<model::actions::TimeStepDependentAction> integrator) {
        scheme->_integrator = std::move(integrator);
        return *this;
    }

    /**
     * Instruct the configurator that the simulation loop should use the provided integrator type
     * @tparam IntegratorType the integrator type
     * @return reference to self
     */
    SchemeConfigurator &withEulerBDIntegrator() {
        scheme->_integrator = scheme->kernel->actions().eulerBDIntegrator(0.);
        return *this;
    }

    /**
     * Instruct the configurator that the simulation loop should use the provided integrator type
     * @param integratorName the integrator type
     * @return reference to self
     */
    SchemeConfigurator &withIntegrator(const std::string &integratorName) {
        scheme->_integrator = scheme->kernel->actions().createIntegrator(integratorName, 0.);
        return *this;
    }

    /**
     * Instruct the configurator that the simulation loop should use the provided reaction scheduler
     * @tparam ReactionSchedulerType the reaction scheduler type
     * @return reference to self
     */
    template<typename ReactionSchedulerType>
    SchemeConfigurator &withReactionScheduler() {
        return withReactionScheduler(detail::identity<ReactionSchedulerType>());
    }

    /**
     * Instruct the configurator that the simulation loop should use the provided reaction scheduler instance.
     * @param scheduler the reaction scheduler instance
     * @return reference to self
     */
    SchemeConfigurator &withReactionScheduler(std::shared_ptr<model::actions::TimeStepDependentAction> scheduler) {
        scheme->_reactions = std::move(scheduler);
        return *this;
    }

    /**
     * Instruct the configurator that the simulation loop should use the provided reaction scheduler
     * @param name the reaction scheduler type
     * @return reference to self
     */
    SchemeConfigurator &withReactionScheduler(const std::string &name) {
        scheme->_reactions = scheme->kernel->actions().createReactionScheduler(name, 0.);
        return *this;
    }

    /**
     * Whether topology reactions should be evaluated.
     * @param evaluate true if they should be evaluated
     * @return reference to self
     */
    SchemeConfigurator &evaluateTopologyReactions(bool evaluate = true) {
        if(evaluate) {
            scheme->_topologyReactions = scheme->kernel->actions().evaluateTopologyReactions(c_::zero);
        } else {
            scheme->_topologyReactions = nullptr;
        }
        return *this;
    }

    /**
     * Whether observables should be evaluated. Also affects wheter trajectories are written.
     * @param evaluate true if they should be evaluated
     * @return reference to self
     */
    SchemeConfigurator &evaluateObservables(bool evaluate = true) {
        scheme->_evaluateObservables = evaluate;
        evaluateObservablesSet = true;
        return *this;
    }

    /**
     * Whether forces should be calculated
     * @param include true if they should be calculated
     * @return reference to self
     */
    SchemeConfigurator &includeForces(bool include = true) {
        if (include) {
            scheme->_forces = scheme->kernel->actions().calculateForces();
        } else {
            scheme->_forces = nullptr;
        }
        includeForcesSet = true;
        return *this;
    }

    /**
     * Instructs to write the system's configuration into the provided file upon simulation start.
     * @param file the file
     * @return reference to self
     */
    SchemeConfigurator &writeConfigToFile(File& file) {
        scheme->configGroup = std::make_unique<h5rd::Group>(file.createGroup("readdy/config"));
        return *this;
    }

    /**
     * Parameter for the neighbor list's skin size, in particular this also increases the box size in cell-linked lists
     * and therefore can be used to reduce memory requirements.
     * @param skin the skin size
     * @return reference to self
     */
    SchemeConfigurator &withSkinSize(scalar skin = 0) {
        skinSize = skin;
        return *this;
    }

    /**
     * Instantiates the scheme with the configuration.
     * @param timeStep the time step
     * @return a unique pointer to the scheme instance
     */
    std::unique_ptr<ReaDDyScheme> configure(scalar timeStep, bool checkTimeStep = true) {

        {
            const double threshold = .1;
            const auto &reactionRegistry = scheme->kernel->context().reactions();
            if(checkTimeStep) {
                for (const auto &o1Reaction : reactionRegistry.order1Flat()) {
                    if (o1Reaction->rate() * timeStep > threshold) {
                        log::warn("Specified a reaction with rate {} and a timestep of {}. Thus, the inverse "
                                  "reaction rate is faster than the time step by at least a factor of {}.",
                                  o1Reaction->rate(), timeStep, 1. / threshold);
                    }
                }
                for (const auto &o2Reaction : reactionRegistry.order2Flat()) {
                    if (o2Reaction->rate() * timeStep > threshold) {
                        log::warn("Specified a reaction with rate {} and a timestep of {}. Thus, the inverse "
                                  "reaction rate is faster than the time step by at least a factor of {}.",
                                  o2Reaction->rate(), timeStep, 1. / threshold);
                    }
                }
            }
        }

        using default_integrator = readdy::model::actions::EulerBDIntegrator;
        using default_reactions = readdy::model::actions::reactions::Gillespie;
        using calculate_forces = readdy::model::actions::CalculateForces;
        using update_neighbor_list = readdy::model::actions::UpdateNeighborList;
        if (useDefaults) {
            if (!scheme->_integrator) {
                scheme->_integrator = scheme->kernel->actions().eulerBDIntegrator(timeStep);
            }
            if (!scheme->_reactions) {
                scheme->_reactions = scheme->kernel->actions().gillespie(timeStep);
            }
            if (!evaluateObservablesSet) {
                scheme->_evaluateObservables = true;
            }
            if (!scheme->_forces && !includeForcesSet) {
                scheme->_forces = scheme->kernel->actions().calculateForces();
            }
        }
        if (scheme->_forces || scheme->_reactions) {
            using ops = model::actions::UpdateNeighborList::Operation;
            scheme->_initNeighborList = scheme->kernel->actions().updateNeighborList(ops::init, skinSize);
            scheme->_neighborList = scheme->kernel->actions().updateNeighborList(ops::update, skinSize);
            scheme->_clearNeighborList = scheme->kernel->actions().updateNeighborList(ops::clear, skinSize);
        }
        if (scheme->_integrator) scheme->_integrator->setTimeStep(timeStep);
        if (scheme->_reactions) scheme->_reactions->setTimeStep(timeStep);
        if (scheme->_topologyReactions) scheme->_topologyReactions->setTimeStep(timeStep);
        std::unique_ptr<ReaDDyScheme> ptr = std::move(scheme);
        scheme = nullptr;
        return ptr;
    }

    /**
     * Instantiates and runs the scheme with the configuration.
     * @param steps number of steps
     * @param timeStep the time step
     * @see configure(timeStep)
     */
    void configureAndRun(const time_step_type steps, scalar timeStep, bool checkTimeStep=true) {
        configure(timeStep, checkTimeStep)->run(steps);
    }

protected:

    SchemeConfigurator &withReactionScheduler(detail::identity<model::actions::reactions::Gillespie>) {
        scheme->_reactions = scheme->kernel->actions().gillespie(c_::zero);
        return *this;
    }

    SchemeConfigurator &withReactionScheduler(detail::identity<model::actions::reactions::UncontrolledApproximation>) {
        scheme->_reactions = scheme->kernel->actions().uncontrolledApproximation(c_::zero);
        return *this;
    }

    /**
     * whether to use defaults or start on a blank slate
     */
    bool useDefaults;
    /**
     * whether evaluateObservables was manually set
     */
    bool evaluateObservablesSet = false;
    /**
     * whether includeForces was manually set
     */
    bool includeForcesSet = false;
    /**
     * the resulting scheme instance
     */
    std::unique_ptr<ReaDDyScheme> scheme = nullptr;
    /**
     * the skin size, has to be saved until the neighbor list actions are created
     */
    scalar skinSize = 0;
};

NAMESPACE_END(api)
NAMESPACE_END(readdy)
