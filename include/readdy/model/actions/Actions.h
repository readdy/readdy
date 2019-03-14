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

    ~AddParticles() override = default;

    void perform() override;

protected:
    std::vector<Particle> particles;
    Kernel *const kernel;
};

class EulerBDIntegrator : public TimeStepDependentAction {
public:
    explicit EulerBDIntegrator(scalar timeStep);

    ~EulerBDIntegrator() override = default;
};

class MdgfrdIntegrator : public TimeStepDependentAction {
public:
    explicit MdgfrdIntegrator(scalar timeStep);

    ~MdgfrdIntegrator() override = default;
};

/**
 * Calculates all forces and energies resulting from potentials (external, pair-potentials, bonded).
 * Optionally calculate the virial which is required if the pressure in the system shall be measured.
 */
class CalculateForces : public Action {
public:
    CalculateForces();

    ~CalculateForces() override = default;
};

class CreateNeighborList : public Action {
public:
    explicit CreateNeighborList(scalar cutoffDistance);

    ~CreateNeighborList() override = default;

    scalar &cutoffDistance() { return _cutoffDistance; }

    const scalar &cutoffDistance() const { return _cutoffDistance; }

protected:
    scalar _cutoffDistance;
};

class UpdateNeighborList : public Action {};

class ClearNeighborList : public Action {};

NAMESPACE_BEGIN(reactions)

class UncontrolledApproximation : public TimeStepDependentAction {

public:
    explicit UncontrolledApproximation(scalar timeStep);

    ~UncontrolledApproximation() override = default;
};

class Gillespie : public TimeStepDependentAction {
public:
    explicit Gillespie(scalar timeStep);

    ~Gillespie() override = default;
};


class DetailedBalance : public TimeStepDependentAction {
public:
    explicit DetailedBalance(scalar timeStep);

    ~DetailedBalance() override = default;
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

    ~EvaluateTopologyReactions() override = default;
};

NAMESPACE_END(top)

/**
 * The Compartments feature defines compartments via characteristic functions that map from Vec3 to bool.
 * For every compartment one can then define conversions that should take place as soon as a particle
 * enters the compartment. Note that the user is responsible for keeping the compartments disjoint.
 *
 * The EvaluateCompartments action performs these conversions.
 */
class EvaluateCompartments : public Action {
public:
    explicit EvaluateCompartments() : Action() {}

    ~EvaluateCompartments() override = default;
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
const std::string getActionName(typename std::enable_if<std::is_base_of<MdgfrdIntegrator, T>::value>::type * = 0) {
    return "MdgfrdIntegrator";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<CalculateForces, T>::value>::type * = 0) {
    return "Calculate forces";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<CreateNeighborList, T>::value>::type * = 0) {
    return "Create neighbor list";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<UpdateNeighborList, T>::value>::type * = 0) {
    return "Update neighbor list";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<ClearNeighborList, T>::value>::type * = 0) {
    return "Clear neighbor list";
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
