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

#if READDY_OSX || READDY_WINDOWS
#include <functional>
#endif

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(actions)

class AddParticles : public Action {

public:
    AddParticles(Kernel *kernel, const std::vector<Particle> &particles);
    AddParticles(Kernel *kernel, const Particle& particle);

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

    virtual bool supportsSkin() const = 0;


protected:
    const Operation operation;
    const scalar skinSize;
};

NAMESPACE_BEGIN(reactions)

class UncontrolledApproximation : public TimeStepDependentAction {

public:
    using reaction_11 = std::function<model::Particle(const model::Particle &)>;
    using reaction_12 = std::function<void(const model::Particle &, model::Particle &, model::Particle &)>;
    using reaction_21 = std::function<model::Particle(const model::Particle &, const model::Particle &)>;
    using reaction_22 = std::function<void(const model::Particle &, const model::Particle &, model::Particle &,
                                           model::Particle &)>;

    explicit UncontrolledApproximation(scalar timeStep);

    virtual void registerReactionScheme_11(const std::string &reactionName, reaction_11 fun) = 0;

    virtual void registerReactionScheme_12(const std::string &reactionName, reaction_12 fun) = 0;

    virtual void registerReactionScheme_21(const std::string &reactionName, reaction_21 fun) = 0;

    virtual void registerReactionScheme_22(const std::string &reactionName, reaction_22 fun) = 0;

};

class Gillespie : public TimeStepDependentAction {
public:
    explicit Gillespie(scalar timeStep);
};

/**
 * Types of reversible reactions
 * - FusionFission A + B <--> C
 * - ConversionConversion A <--> B
 * - EnzymaticEnzymatic A + C <--> B + C, reaction radii forw. and backw. have to be equal!
 */
enum ReversibleType {
    FusionFission, ConversionConversion, EnzymaticEnzymatic
};

inline std::ostream& operator<<(std::ostream& os, const ReversibleType& reversibleType) {
    switch (reversibleType) {
        case ReversibleType::FusionFission: os << "FusionFission"; break;
        case ReversibleType::ConversionConversion: os << "ConversionConversion"; break;
        case ReversibleType::EnzymaticEnzymatic: os << "EnzymaticEnzymatic"; break;
    }
    return os;
}

/**
 * A reversible reaction `lhs <--> rhs` is defined by two reactions:
 * - forward `lhs --> rhs`
 * - backward `lhs <-- rhs`
 */
struct ReversibleReactionConfig {
    using ReactionId = readdy::model::reactions::Reaction::reaction_id;
    using pot = readdy::model::potentials::PotentialOrder2 *;

    ReactionId forwardId;
    ReactionId backwardId;
    const model::reactions::Reaction *forwardReaction;
    const model::reactions::Reaction *backwardReaction;
    std::string forwardName;
    std::string backwardName;

    std::uint8_t numberLhsTypes; // either 1 or 2
    std::array<particle_type_type, 2> lhsTypes;
    std::array<std::string, 2> lhsNames;
    std::uint8_t numberRhsTypes; // either 1 or 2
    std::array<particle_type_type, 2> rhsTypes;
    std::array<std::string, 2> rhsNames;

    ReversibleType reversibleType;

    scalar microForwardRate;
    scalar microBackwardRate;

    scalar reactionRadius = 0.; // must be equal forward and backward only FusionFission and EnzymaticEnzymatic

    std::vector<pot> lhsPotentials;
    std::vector<pot> rhsPotentials; // only needed for EnzymaticEnzymatic
    scalar lhsInteractionRadius;
    scalar rhsInteractionRadius;
    scalar lhsInteractionVolume;
    scalar rhsInteractionVolume;
    scalar effectiveLhsInteractionVolume;
    scalar effectiveRhsInteractionVolume;
    scalar effectiveLhsReactionVolume;
    scalar effectiveRhsReactionVolume;

    scalar totalVolume;
    scalar kbt;
    scalar acceptancePrefactor = 1.; // EnzymaticEnzymatic only

    // To draw fission radius
    std::vector<scalar> fissionRadii; // FusionFission only
    std::vector<scalar> cumulativeFissionProb; // FusionFission only

    // these are inferred from the microscopic quantities
    scalar equilibriumConstant;
    scalar macroForwardRate;
    scalar macroBackwardRate;

    explicit ReversibleReactionConfig(ReactionId forwardId, ReactionId backwardId, const Context &ctx);

    std::string describe() const;

    scalar drawFissionDistance() const {
        auto u = readdy::model::rnd::uniform_real();
        auto it = std::lower_bound(cumulativeFissionProb.begin(), cumulativeFissionProb.end(), u);
        auto index = std::distance(cumulativeFissionProb.begin(), it);
        return fissionRadii[index];
    }
};

/**
 * Comparator for identifying equivalent reversible reactions, which are forbidden due to ambiguity.
 */
inline const bool equivalentReversibleReactions(const ReversibleReactionConfig &rev1, const ReversibleReactionConfig &rev2) {
    return (rev1.lhsTypes == rev2.lhsTypes and rev1.rhsTypes == rev2.rhsTypes) or
           (rev1.lhsTypes == rev2.rhsTypes and rev1.rhsTypes == rev2.lhsTypes);
};

class DetailedBalance : public TimeStepDependentAction {
public:
    explicit DetailedBalance(scalar timeStep);

    const std::vector<std::shared_ptr<const ReversibleReactionConfig>> &reversibleReactions() const {
        return _reversibleReactionsContainer;
    }
    
    std::string describe() const;

protected:
    void searchReversibleReactions(const Context& ctx);

    std::vector<std::shared_ptr<const ReversibleReactionConfig>> _reversibleReactionsContainer;
    // the map provides a view on the container, where two (unidirectional) reaction ids
    // might point to the same ReversibleReactionConfig
    std::unordered_map<model::reactions::Reaction::reaction_id, std::shared_ptr<const ReversibleReactionConfig>>
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
const std::string getActionName(typename std::enable_if<std::is_base_of<reactions::DetailedBalance, T>::value>::type * = 0) {
    return "DetailedBalance";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<top::EvaluateTopologyReactions, T>::value>::type * = 0) {
    return "EvaluateTopologyReactions";
}

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<EvaluateCompartments, T>::value>::type * = 0) {
    return "EvaluateCompartments";
}

NAMESPACE_END(actions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
