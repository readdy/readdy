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
 * This files contains a selection of possible Actions that can be executed on a kernel:
 *   - AddParticles: An action with which particles can be added.
 *   - EulerBDIntegrator: Propagates the particles through the system.
 *   - UpdateNeighborList: Creates and updates neighbor lists.
 *   - CalculateForces: Calculates forces for later use in, e.g., integration schemes.
 *   - UncontrolledApproximation: Executes reactions, resolving conflicts in an uncontrolled way.
 *   - Gillespie: Executes reactions, sampling one reaction event after the other weighted with their rates.
 *   - DetailedBalance: Executes reactions, and assures detailed balance for reversible reactions.
 *   - EvaluateCompartments: Perform instantaneous particle conversions depending on the particles' position.
 *   - EvaluateTopologyReactions: Execute reactions involving topologies.
 *
 * Further, specializations of ActionName<T> are declared.
 *
 * @file Actions.h
 * @brief Declaration of all globally available actions.
 * @author clonker
 * @author chrisfroe
 * @date 11.04.16
 * @todo provide more detailed descriptions for some of the actions
 */
#pragma once

#include <type_traits>
#include <readdy/model/Particle.h>
#include <readdy/model/actions/Action.h>
#include <readdy/model/observables/Observable.h>
#include <readdy/model/reactions/Reaction.h>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <readdy/model/Context.h>
#include <readdy/model/actions/DetailedBalance.h>

#if READDY_OSX || READDY_WINDOWS
#include <functional>
#endif

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(actions)

class AddParticles : public Action {

public:
    AddParticles(Kernel *kernel, const std::vector<Particle> &particles);

    void perform(const util::PerformanceNode &node) override;

protected:
    std::vector<Particle> particles;
    Kernel *const kernel;
};

class EulerBDIntegrator : public TimeStepDependentAction {
public:
    explicit EulerBDIntegrator(scalar timeStep);
};

class CalculateForces : public Action {
public:
    CalculateForces();
};

class UpdateNeighborList : public Action {
public:
    enum Operation {
        init, update, clear
    };

    explicit UpdateNeighborList(Operation operation = Operation::init, scalar skinSize = 0);

    scalar &skin() { return skinSize; }

    const scalar &skin() const { return skinSize; }

protected:
    Operation operation;
    scalar skinSize;
};

NAMESPACE_BEGIN(reactions)

class UncontrolledApproximation : public TimeStepDependentAction {

public:
    explicit UncontrolledApproximation(scalar timeStep);
};

class Gillespie : public TimeStepDependentAction {
public:
    explicit Gillespie(scalar timeStep);
};


class DetailedBalance : public TimeStepDependentAction {
public:
    explicit DetailedBalance(scalar timeStep);

    const std::vector<std::shared_ptr<const ReversibleReactionConfig>> &reversibleReactions() const {
        return _reversibleReactionsContainer;
    }

    std::string describe() const;

protected:
    void searchReversibleReactions(const Context &ctx);

    std::vector<std::shared_ptr<const ReversibleReactionConfig>> _reversibleReactionsContainer;
    // the map provides a view on the container for quick runtime lookup,
    // usually two (unidirectional) reaction ids point to the same ReversibleReactionConfig
    std::unordered_map<model::reactions::Reaction::ReactionId, std::shared_ptr<const ReversibleReactionConfig>>
            _reversibleReactionsMap;
};

NAMESPACE_END(reactions)

NAMESPACE_BEGIN(top)

class EvaluateTopologyReactions : public TimeStepDependentAction {
public:
    explicit EvaluateTopologyReactions(scalar timeStep);
};

NAMESPACE_END(top)

class EvaluateCompartments : public Action {
public:
    explicit EvaluateCompartments() : Action() {}
};

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<AddParticles, T>::value>::type * = 0) {
    return "AddParticles";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<EulerBDIntegrator, T>::value>::type * = 0) {
    return "EulerBDIntegrator";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<CalculateForces, T>::value>::type * = 0) {
    return "Calculate forces";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<UpdateNeighborList, T>::value>::type * = 0) {
    return "Update neighbor list";
}

template<typename T>
const std::string
getActionName(typename std::enable_if<std::is_base_of<reactions::UncontrolledApproximation, T>::value>::type * = 0) {
    return "UncontrolledApproximation";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<reactions::Gillespie, T>::value>::type * = 0) {
    return "Gillespie";
}

template<typename T>
const std::string
getActionName(typename std::enable_if<std::is_base_of<reactions::DetailedBalance, T>::value>::type * = 0) {
    return "DetailedBalance";
}

template<typename T>
const std::string
getActionName(typename std::enable_if<std::is_base_of<top::EvaluateTopologyReactions, T>::value>::type * = 0) {
    return "EvaluateTopologyReactions";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<EvaluateCompartments, T>::value>::type * = 0) {
    return "EvaluateCompartments";
}

NAMESPACE_END(actions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
