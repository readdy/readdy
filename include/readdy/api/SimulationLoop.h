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
 * Vital parts of the loop api are defined in this header. First of all, there is the SimulationScheme. Its task is
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
 * @file SimulationLoop.h
 * @brief Header file containing the SchemeConfigurator<T> and a corresponding SimulationLoop superclass definition.
 * @author clonker
 * @author chrisfroe
 * @date 23.08.16
 */

#pragma once

#include <memory>
#include <utility>
#include <type_traits>

#include <h5rd/h5rd.h>

#include <readdy/common/common.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/IOUtils.h>

#include "Saver.h"

namespace readdy::api {

/**
 * superclass for all simulation loops
 */
class SimulationLoop {
public:
    /**
     * the type of function that is responsible for deciding whether the simulation should continue
     */
    using continue_fun = std::function<bool(TimeStep)>;

    using TimeStepActionPtr = std::shared_ptr<readdy::model::actions::TimeStepDependentAction>;
    using ActionPtr = std::shared_ptr<readdy::model::actions::Action>;

    /**
     * Creates a new simulation scheme. Creates and initializes actions: Sets the neighborlist distance
     * to be the largest cutoff present in the context.
     * @param kernel the kernel
     * @param timeStep the time step width
     */
    explicit SimulationLoop(model::Kernel *const kernel, scalar timeStep)
            : _kernel(kernel), _timeStep(timeStep),
              _integrator(kernel->actions().eulerBDIntegrator(timeStep).release()),
              _reactions(
                      kernel->supportsGillespie() ?
                      static_cast<readdy::model::actions::TimeStepDependentAction *>(
                              kernel->actions().gillespie(timeStep).release()) :
                      static_cast<readdy::model::actions::TimeStepDependentAction *>(
                              kernel->actions().uncontrolledApproximation(timeStep).release())),
              _forces(kernel->actions().calculateForces().release()),
              _initNeighborList(kernel->actions().createNeighborList(kernel->context().calculateMaxCutoff()).release()),
              _updateNeighborList(kernel->actions().updateNeighborList().release()),
              _clearNeighborList(kernel->actions().clearNeighborList().release()),
              _topologyReactions(kernel->supportsTopologies() ?
                                 kernel->actions().evaluateTopologyReactions(timeStep).release() : nullptr) {}

    /**
     * This function gives access to an progressCallback function that gets called every 100 time steps if one is
     * running the simulation for a fixed number of time steps.
     * @return the callback function
     */
    std::function<void(TimeStep)> &progressCallback() {
        return _progressCallback;
    }

    const std::function<void(TimeStep)> &progressCallback() const {
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

    std::size_t checkpointingStride() const { return _checkpointingStride; }
    void setCheckpointingStride(std::size_t stride) { _checkpointingStride = stride; }


    void runInitialize() {
        _kernel->initialize();
        if (configGroup) {
            model::ioutils::writeSimulationSetup(*configGroup, _kernel->context());
        }
    }

    void runInitializeNeighborList() {
        if (_initNeighborList) _initNeighborList->perform();
    }

    void runUpdateNeighborList() {
        if (_updateNeighborList) _updateNeighborList->perform();
    }

    void runClearNeighborList() {
        if (_clearNeighborList) _clearNeighborList->perform();
    }

    void runForces() {
        if (_forces) _forces->perform();
    }

    void runEvaluateObservables(TimeStep t) {
        if (_evaluateObservables) _kernel->evaluateObservables(t);
    }

    void runIntegrator() {
        if (_integrator) _integrator->perform();
    }

    void runReactions() {
        if (_reactions) _reactions->perform();
    }

    void runTopologyReactions() {
      if (_topologyReactions) {
          _topologyReactions->perform();
      }
    }

    TimeStepActionPtr &integrator() { return _integrator; }

    const TimeStepActionPtr &integrator() const { return _integrator; }

    void useIntegrator(const std::string &name, scalar timeStep = -1) {
        _integrator = _kernel->actions().createIntegrator(name, timeStep > 0 ? timeStep : _timeStep);
    }

    TimeStepActionPtr &reactionScheduler() { return _reactions; }

    const TimeStepActionPtr &reactionScheduler() const { return _reactions; }

    void useReactionScheduler(const std::string &name, scalar timeStep = -1) {
        _reactions = _kernel->actions().createReactionScheduler(name, timeStep > 0 ? timeStep : _timeStep);
    }

    void evaluateTopologyReactions(bool evaluate, scalar timeStep = -1) {
        _topologyReactions = evaluate ? _kernel->actions().evaluateTopologyReactions(
                timeStep > 0 ? timeStep : _timeStep) : nullptr;
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

    scalar &neighborListCutoff() {
        return _initNeighborList->cutoffDistance();
    }

    [[nodiscard]] scalar neighborListCutoff() const {
        return _initNeighborList->cutoffDistance();
    }

    /**
     * This method runs a simulation for a fixed number of time steps by providing an appropriate continue function.
     * @param steps the number of steps
     */
    virtual void run(TimeStep steps) {
        // show every 100 time steps
        if (!_progressCallback) {
            _progressCallback = [this, steps](TimeStep current) {
                log::info("Simulation progress: {} / {} steps", (current - _start), steps);
            };
        }
        auto defaultContinueCriterion = [this, steps](const TimeStep current) {
            if (current != _start && _progressOutputStride > 0 && (current - _start) % _progressOutputStride == 0) {
                _progressCallback(current);
            }
            return current < _start + steps;
        };
        run(defaultContinueCriterion);
    }

    /**
     * ReaDDy scheme implementation of the simulation loop
     * @param continueFun the continue function
     */
    void run(const continue_fun &continueFun) {
        validate(_timeStep);
        {
            bool requiresNeighborList = _initNeighborList->cutoffDistance() > 0;
            if (requiresNeighborList) {
                if (!(_initNeighborList && _updateNeighborList && _clearNeighborList)) {
                    throw std::logic_error("Neighbor list required but set to null!");
                }
                _initNeighborList->cutoffDistance() = std::max(_initNeighborList->cutoffDistance(),
                                                               kernel()->context().calculateMaxCutoff());
            }
            runInitialize();
            if (requiresNeighborList) runInitializeNeighborList();
            runForces();
            TimeStep t = _start;
            if(_saver) {
                // this needs to happen before observables because observables can in principle influence the state
                _saver->makeCheckpoint(_kernel, t);
            }
            runEvaluateObservables(t);
            std::for_each(std::begin(_callbacks), std::end(_callbacks), [t](const auto &callback) {
                callback(t);
            });
            while (continueFun(t)) {
                runIntegrator();
                if (requiresNeighborList) runUpdateNeighborList();
                runReactions();
                runTopologyReactions();
                if (requiresNeighborList) runUpdateNeighborList();
                runForces();
                if(_saver && (t + 1) % _checkpointingStride == 0) {
                    // this needs to happen before observables because observables can in principle influence the state
                    _saver->makeCheckpoint(_kernel, t + 1);
                }
                runEvaluateObservables(t + 1);
                std::for_each(std::begin(_callbacks), std::end(_callbacks), [t](const auto &callback) {
                    callback(t + 1);
                });
                ++t;

                _kernel->stateModel().setTime(_kernel->stateModel().time() + _timeStep);
            }
            if (requiresNeighborList) runClearNeighborList();
            _start = t;
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

    void addCallback(std::function<void(TimeStep)> f) {
        _callbacks.emplace_back(std::move(f));
    }

    void setSaver(std::shared_ptr<Saver> saver) {
        _saver = std::move(saver);
    }

    std::shared_ptr<Saver> saver() const {
        return _saver;
    }

    std::string describe() {
        namespace rus = readdy::util::str;
        std::string description;
        description += fmt::format("Configured simulation loop with:\n");
        description += fmt::format("--------------------------------\n");
        description += fmt::format(" - timeStep = {}\n", _timeStep);
        description += fmt::format(" - evaluateObservables = {}\n", _evaluateObservables);
        description += fmt::format(" - progressOutputStride = {}\n", _progressOutputStride);
        description += fmt::format(" - context written to file = {}\n", static_cast<bool>(configGroup));
        // todo let actions know their name?
        description += fmt::format(" - Performing actions:\n");
        description += fmt::format("   * Initialize neighbor list? {}\n", static_cast<bool>(_initNeighborList));
        description += fmt::format("   * Update neighbor list? {}\n", static_cast<bool>(_updateNeighborList));
        description += fmt::format("   * Clear neighbor list? {}\n", static_cast<bool>(_clearNeighborList));
        description += fmt::format("   * Integrate diffusion? {}\n", static_cast<bool>(_integrator));
        description += fmt::format("   * Calculate forces? {}\n", static_cast<bool>(_forces));
        description += fmt::format("   * Handle reactions? {}\n", static_cast<bool>(_reactions));
        description += fmt::format("   * Handle topology reactions? {}\n", static_cast<bool>(_topologyReactions));
        if (_saver) {
            description += fmt::format(" - Performing checkpointing:\n");
            description += fmt::format("   * stride: {}\n", _checkpointingStride);
            description += fmt::format("   * base path: {}\n", _saver->basePath());
            description += fmt::format("   * checkpoint filename template: {}\n", _saver->checkpointTemplate());
            description += fmt::format("   * maximal number saves: {}\n", _saver->maxNSaves());
        }
        return description;
    }

protected:
    model::Kernel *const _kernel;
    std::shared_ptr<model::actions::TimeStepDependentAction> _integrator{nullptr};
    std::shared_ptr<model::actions::Action> _forces{nullptr};
    std::shared_ptr<model::actions::TimeStepDependentAction> _reactions{nullptr};
    std::shared_ptr<model::actions::CreateNeighborList> _initNeighborList{nullptr};
    std::shared_ptr<model::actions::UpdateNeighborList> _updateNeighborList{nullptr};
    std::shared_ptr<model::actions::top::EvaluateTopologyReactions> _topologyReactions{nullptr};
    std::shared_ptr<model::actions::ClearNeighborList> _clearNeighborList{nullptr};
    std::shared_ptr<api::Saver> _saver {nullptr};
    std::shared_ptr<h5rd::Group> configGroup{nullptr};

    bool _evaluateObservables = true;
    TimeStep _start = 0;
    std::size_t _progressOutputStride = 100;
    std::size_t _checkpointingStride = 10000;
    std::function<void(TimeStep)> _progressCallback;
    scalar _timeStep;

    std::vector<std::function<void(TimeStep)>> _callbacks;
};

}
