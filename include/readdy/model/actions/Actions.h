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
 *   - CreateNeighborList: Creates neighbor lists.
 *   - UpdateNeighborList: Updates neighbor lists.
 *   - ClearNeighborList: Clears neighbor lists.
 *   - CalculateForces: Calculates forces for later use in, e.g., integration schemes.
 *   - UncontrolledApproximation: Executes reactions, resolving conflicts in an uncontrolled way.
 *   - Gillespie: Executes reactions, sampling one reaction event after the other weighted with their rates.
 *   - DetailedBalance: Executes reactions, and assures detailed balance for reversible reactions.
 *   - EvaluateCompartments: Perform instantaneous particle conversions depending on the particles' position.
 *   - EvaluateTopologyReactions: Execute reactions involving topologies.
 *   - BreakBonds: Remove edges in topologies based on current bond-energy.
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
#include <readdy/common/index_persistent_vector.h>
#include "Utils.h"
#include <readdy/model/topologies/GraphTopology.h>

#if READDY_OSX || READDY_WINDOWS
#include <functional>
#endif

namespace readdy::model::actions {

class AddParticles : public Action {

public:
    AddParticles(Kernel *kernel, std::vector<Particle> particles);

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

class UpdateNeighborList : public Action {
};

class ClearNeighborList : public Action {
};

namespace reactions {

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
    std::unordered_map<ReactionId, std::shared_ptr<const ReversibleReactionConfig>> _reversibleReactionsMap;
};

}

namespace top {

class EvaluateTopologyReactions : public TimeStepDependentAction {
public:
    explicit EvaluateTopologyReactions(scalar timeStep);

    ~EvaluateTopologyReactions() override = default;
};

class BreakConfig {
public:
    void addBreakablePair(ParticleTypeId type1, ParticleTypeId type2, scalar thresholdEnergy, scalar rate) {
        if (rate <= 0.) {
            throw std::logic_error(fmt::format("Rate of breaking a bond must be positive, was {}", rate));
        }
        _thresholdEnergies.emplace(std::make_pair(std::make_pair(type1, type2), thresholdEnergy));
        _breakRates.emplace(std::make_pair(std::make_pair(type1, type2), rate));
    }

    const util::particle_type_pair_unordered_map<scalar> &thresholdEnergies() const {
        return _thresholdEnergies;
    }

    const util::particle_type_pair_unordered_map<scalar> &breakRates() const {
        return _breakRates;
    }

private:
    util::particle_type_pair_unordered_map<scalar> _thresholdEnergies{};
    util::particle_type_pair_unordered_map<scalar> _breakRates{};
};

/**
 * BreakBonds can remove edges in topology graphs, based on the current potential energy that the bond
 * (belonging to the edge) contains. Uses the breakConfig, which holds threshold energies and rates (frequencies)
 * at which edge-removals are attempted. These values are defined per pair of particle types.
 */
class BreakBonds : public TimeStepDependentAction {
public:
    using vertex_ref = readdy::model::top::Graph::VertexIndex;

    explicit BreakBonds(scalar timeStep, BreakConfig breakConfig);

    ~BreakBonds() override = default;

    const util::particle_type_pair_unordered_map<scalar> &thresholdEnergies() const {
        return breakConfig.thresholdEnergies();
    }

    const util::particle_type_pair_unordered_map<scalar> &breakRates() const {
        return breakConfig.breakRates();
    }

protected:
    const BreakConfig breakConfig;

    template<typename Kernel, typename TopologyRef, typename Model, typename ParticleData>
    void genericPerform(readdy::util::index_persistent_vector<TopologyRef> &topologies, Model &model, Kernel *kernel, ParticleData &particleData) {
        std::vector<readdy::model::top::GraphTopology> resultingTopologies;
        std::size_t topologyIdx = 0;
        for (auto &top : topologies) {
            if (!top->isDeactivated()) {
                auto reactionFunction = [&](
                        readdy::model::top::GraphTopology &t) -> readdy::model::top::reactions::Recipe {
                    readdy::model::top::reactions::Recipe recipe(t);
                    for (const auto &edge : t.graph().edges()) {
                        auto [e1, e2] = edge;
                        auto energy = evaluateEdgeEnergy(edge, t, kernel);
                        auto v1Type = t.graph().vertices().at(e1)->particleType;
                        auto v2Type = t.graph().vertices().at(e2)->particleType;
                        const auto thresholdEnergyIt = thresholdEnergies().find(std::tie(v1Type, v2Type));
                        if (thresholdEnergyIt != thresholdEnergies().end()) {
                            if (energy > thresholdEnergyIt->second) {
                                const auto &rate = breakRates().at(std::tie(v1Type, v2Type));
                                if (readdy::model::rnd::uniform_real() < 1 - std::exp(-rate * _timeStep)) {
                                    recipe.removeEdge(edge);
                                }
                            }
                        }
                    }
                    return std::move(recipe);
                };
                scalar rateDoesntMatter{1.};
                readdy::model::top::reactions::StructuralTopologyReaction reaction("__internal_break_bonds",
                                                                                   reactionFunction, rateDoesntMatter);
                readdy::model::actions::top::executeStructuralReaction(topologies, resultingTopologies, top, reaction,
                                                                       topologyIdx, particleData, kernel);

            }
            ++topologyIdx;
        }

        const auto &context = kernel->context();
        for (auto &&newTopology : resultingTopologies) {
            if (!newTopology.isNormalParticle(*kernel)) {
                // we have a new topology here, update data accordingly.
                newTopology.updateReactionRates(
                        context.topologyRegistry().structuralReactionsOf(newTopology.type()));
                newTopology.configure();
                model.insert_topology(std::move(newTopology));
            } else {
                // if we have a single particle that is not of flavor topology, remove from topology structure!
                auto particleIndex = 0;
                {
                    // find first particle that is not deactivated
                    for(auto it = newTopology.graph().vertices().begin(); it != newTopology.graph().vertices().end(); ++it) {
                        if(!it->deactivated()) {
                            particleIndex = std::distance(newTopology.graph().vertices().begin(), it);
                            break;
                        }
                    }
                }
                model.getParticleData()->entry_at(particleIndex).topology_index = -1;
            }
        }
    }

    template <typename Kernel>
    scalar
    evaluateEdgeEnergy(model::top::Graph::Edge edge, const readdy::model::top::GraphTopology &t, Kernel *kernel) const {
        auto [i1, i2] = edge;
        const auto& v1 = t.graph().vertices().at(i1);
        const auto& v2 = t.graph().vertices().at(i2);

        // find bond configurations for given edge
        std::unordered_map<api::BondType, std::vector<readdy::model::top::pot::BondConfiguration>, readdy::util::hash::EnumClassHash> bondConfigs;
        {
            const auto &potentialConfiguration = kernel->context().topologyRegistry().potentialConfiguration();
            auto it = potentialConfiguration.pairPotentials.find(std::tie(v1->particleType, v2->particleType));
            if (it != potentialConfiguration.pairPotentials.end()) {
                for (const auto &cfg : it->second) {
                    bondConfigs[cfg.type].emplace_back(v1->particleIndex, v2->particleIndex,
                            cfg.forceConstant, cfg.length);
                }
            } else {
                std::ostringstream ss;
                auto p1 = t.particleForVertex(i1);
                auto p2 = t.particleForVertex(i2);

                throw std::invalid_argument(fmt::format("The edge {} ({}) == {} ({}) has no bond configured!",
                        v1->particleIndex, kernel->context().particleTypes().nameOf(p1.type()), v2->particleIndex,
                        kernel->context().particleTypes().nameOf(p2.type())));
            }
        }

        // transform configurations to potential instances
        std::vector<std::unique_ptr<readdy::model::top::Topology::BondedPotential>> bondedPotentials;
        for (const auto &bond : bondConfigs) {
            switch (bond.first) {
                case api::BondType::HARMONIC: {
                    bondedPotentials.push_back(
                            std::make_unique<readdy::model::top::TopologyActionFactory::harmonic_bond>(bond.second));
                    break;
                };
            }
        }

        // create actions, perform and accumulate energies on edge
        auto taf = kernel->getTopologyActionFactory();
        scalar totalEnergyForEdge{0.};
        for (const auto &bondedPot : bondedPotentials) {
            totalEnergyForEdge += bondedPot->createForceAndEnergyAction(taf)->perform(&t);
        }

        return totalEnergyForEdge;
    }
};

}

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

template<typename T>
const std::string getActionName(typename std::enable_if<std::is_base_of<top::BreakBonds, T>::value>::type * = 0) {
    return "BreakBonds";
}
}
