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
class SimulationLoop {
    using NeighborListOps = model::actions::UpdateNeighborList::Operation;
public:
    /**
     * the type of function that is responsible for deciding whether the simulation should continue
     */
    using continue_fun = std::function<bool(time_step_type)>;

    using TimeStepActionPtr = std::shared_ptr<readdy::model::actions::TimeStepDependentAction>;
    using ActionPtr = std::shared_ptr<readdy::model::actions::Action>;

    /**
     * Creates a new simulation scheme.
     * @param kernel the kernel
     * @param performanceRoot performance stuff
     */
    explicit SimulationLoop(model::Kernel *const kernel, scalar timeStep, util::PerformanceNode &performanceRoot)
            : _kernel(kernel), _performanceRoot(performanceRoot), _timeStep(timeStep),
              _integrator(kernel->actions().eulerBDIntegrator(timeStep).release()),
              _reactions(kernel->actions().gillespie(timeStep).release()),
              _forces(kernel->actions().calculateForces().release()),
              _initNeighborList(kernel->actions().updateNeighborList(NeighborListOps::init).release()),
              _neighborList(kernel->actions().updateNeighborList(NeighborListOps::update).release()),
              _clearNeighborList(kernel->actions().updateNeighborList(NeighborListOps::clear).release()),
              _topologyReactions(kernel->actions().evaluateTopologyReactions(timeStep).release()) {}

    /**
     * This function gives access the an updateCallback function that gets called every 100 time steps if one is
     * running the simulation for a fixed number of time steps.
     * @return the callback function
     */
    std::function<void(time_step_type)> &progressCallback() {
        return _progressCallback;
    }

    const std::function<void(time_step_type)> &progressCallback() const {
        return _progressCallback;
    }

    const std::size_t &progressOutputStride() const {
        return _progressOutputStride;
    }

    /**
     * show progress every N steps
     * @return N steps
     */
    std::size_t &progressOutputStride() {
        return _progressOutputStride;
    }


    void runInitialize() {
        _kernel->initialize();
        if (configGroup) {
            model::ioutils::writeSimulationSetup(*configGroup, _kernel->context());
        }
    }

    void runInitializeNeighborList() {
        if (_initNeighborList) _initNeighborList->perform(_performanceRoot.subnode("initNeighborList"));
    }

    void runUpdateNeighborList() {
        if (_neighborList) _neighborList->perform(_performanceRoot.subnode("neighborList"));
    }

    void runClearNeighborList() {
        if (_clearNeighborList) _clearNeighborList->perform(_performanceRoot.subnode("clearNeighborList"));
    }

    void runForces() {
        if (_forces) _forces->perform(_performanceRoot.subnode("forces"));
    }

    void runEvaluateObservables(time_step_type t) {
        if (_evaluateObservables) _kernel->evaluateObservables(t);
    }

    void runIntegrator() {
        if (_integrator) _integrator->perform(_performanceRoot.subnode("integrator"));
    }

    void runReactions() {
        if (_reactions) _reactions->perform(_performanceRoot.subnode("reactionScheduler"));
    }

    void runTopologyReactions() {
        if (_topologyReactions) _topologyReactions->perform(_performanceRoot.subnode("evaluateTopologyReactions"));
    }

    TimeStepActionPtr &integrator() { return _integrator; }

    const TimeStepActionPtr &integrator() const { return _integrator; }

    void useIntegrator(const std::string &name, scalar timeStep=-1) {
        _integrator = _kernel->actions().createIntegrator(name, timeStep > 0 ? timeStep : _timeStep);
    }

    TimeStepActionPtr &reactionScheduler() { return _reactions; }

    const TimeStepActionPtr &reactionScheduler() const { return _reactions; }

    void useReactionScheduler(const std::string &name, scalar timeStep=-1) {
        _reactions = _kernel->actions().createReactionScheduler(name, timeStep > 0 ? timeStep : _timeStep);
    }

    void evaluateTopologyReactions(bool evaluate, scalar timeStep=-1) {
        _topologyReactions = evaluate ? _kernel->actions().evaluateTopologyReactions(timeStep > 0 ? timeStep : _timeStep) : nullptr;
    }

    void evaluateObservables(bool evaluate) {
        _evaluateObservables = evaluate;
    }

    void evaluateForces(bool include) {
        _forces = include ? _kernel->actions().calculateForces() : nullptr;
    }

    void writeConfigToFile(File &file) {
        configGroup = std::make_unique<h5rd::Group>(file.createGroup("readdy/config"));
    }

    scalar &skinSize() { return _skinSize; }

    const scalar &skinSize() const { return _skinSize; }

    /**
     * This method runs a simulation for a fixed number of time steps by providing an appropriate continue function.
     * @param steps the number of steps
     */
    virtual void run(time_step_type steps) {
        // show every 100 time steps
        if (!_progressCallback) {
            _progressCallback = [this, steps](time_step_type current) {
                log::info("Simulation progress: {} / {} steps", (current - start), steps);
            };
        }
        auto defaultContinueCriterion = [this, steps](const time_step_type current) {
            if (current != start && _progressOutputStride > 0 && (current - start) % _progressOutputStride == 0) {
                _progressCallback(current);
            }
            return current < start + steps;
        };
        run(defaultContinueCriterion);
    };


    /**
     * ReaDDy scheme implementation of the simulation loop
     * @param continueFun the continue function
     */
    void run(const continue_fun &continueFun) {
        validate(_timeStep);
        {
            bool requiresNeighborList = _forces || _reactions || _topologyReactions;
            if (requiresNeighborList) {
                if (!(_initNeighborList && _neighborList && _clearNeighborList)) {
                    throw std::logic_error("Neighbor list required but set to null!");
                }
                _initNeighborList->skin() = _skinSize;
                _neighborList->skin() = _skinSize;
                _clearNeighborList->skin() = _skinSize;
            }
            auto runTimer = _performanceRoot.timeit();
            runInitialize();
            if (requiresNeighborList) runInitializeNeighborList();
            runForces();
            time_step_type t = start;
            runEvaluateObservables(t);
            while (continueFun(t)) {
                runIntegrator();
                if (requiresNeighborList) runUpdateNeighborList();
                runReactions();
                runTopologyReactions();
                if (requiresNeighborList) runUpdateNeighborList();
                runForces();
                runEvaluateObservables(t + 1);
                ++t;
            }
            if (requiresNeighborList) runClearNeighborList();
            start = t;
            log::info("Simulation completed");
        }
    }

    void validate(scalar timeStep) {
        _kernel->context().validate();
        {
            // validate whether reactions are too fast for dt
            const double threshold = .1;
            const auto &reactionRegistry = _kernel->context().reactions();
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

    scalar timeStep() const {
        return _timeStep;
    }
    
    model::Kernel *const kernel() {
        return _kernel;
    }
    
    const model::Kernel *const kernel() const {
        return _kernel;
    }

protected:
    model::Kernel *const _kernel;
    std::shared_ptr<model::actions::TimeStepDependentAction> _integrator{nullptr};
    std::shared_ptr<model::actions::Action> _forces{nullptr};
    std::shared_ptr<model::actions::TimeStepDependentAction> _reactions{nullptr};
    std::shared_ptr<model::actions::UpdateNeighborList> _initNeighborList{nullptr};
    std::shared_ptr<model::actions::UpdateNeighborList> _neighborList{nullptr};
    std::shared_ptr<model::actions::top::EvaluateTopologyReactions> _topologyReactions{nullptr};
    std::shared_ptr<model::actions::UpdateNeighborList> _clearNeighborList{nullptr};
    scalar _skinSize = 0;
    std::shared_ptr<h5rd::Group> configGroup = nullptr;

    bool _evaluateObservables = true;
    time_step_type start = 0;
    std::size_t _progressOutputStride = 100;
    const util::PerformanceNode &_performanceRoot;
    std::function<void(time_step_type)> _progressCallback;
    scalar _timeStep;
};

NAMESPACE_END(api)
NAMESPACE_END(readdy)
