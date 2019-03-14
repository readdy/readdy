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
 * This file contains the declaration of the action factory. Internally, the factory is simply a map of
 * string -> std::function<Action*()>, which then can get called.
 *
 * @file ActionFactory.h
 * @brief Declaration of the action factory.
 * @author clonker
 * @author chrisfroe
 * @date 08.04.16
 */

#pragma once

#include <type_traits>
#include <unordered_map>
#include <vector>
#include <readdy/model/actions/Action.h>
#include <readdy/model/actions/Actions.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(actions)

class ActionFactory {
public:

    virtual std::vector<std::string> getAvailableActions() const {
        return {
                getActionName<AddParticles>(), getActionName<EulerBDIntegrator>(), getActionName<CalculateForces>(),
                getActionName<CreateNeighborList>(), getActionName<UpdateNeighborList>(),
                getActionName<ClearNeighborList>(), getActionName<reactions::UncontrolledApproximation>(),
                getActionName<reactions::Gillespie>(), getActionName<reactions::DetailedBalance>(),
                getActionName<top::EvaluateTopologyReactions>(), getActionName<MdgfrdIntegrator>()
        };
    }

    /*
     * Convenience stuff
     */
    std::unique_ptr<TimeStepDependentAction> createIntegrator(const std::string &name, scalar timeStep) {
        if (name == getActionName<EulerBDIntegrator>()) {
            return std::unique_ptr<TimeStepDependentAction>(eulerBDIntegrator(timeStep));
        }
        throw std::invalid_argument("Requested integrator " + name + " is not available.");
    }

    std::unique_ptr<TimeStepDependentAction> createReactionScheduler(const std::string &name, scalar timeStep) {
        if (name == getActionName<reactions::Gillespie>()) {
            return std::unique_ptr<TimeStepDependentAction>(gillespie(timeStep));
        }
        if (name == getActionName<reactions::UncontrolledApproximation>()) {
            return std::unique_ptr<TimeStepDependentAction>(uncontrolledApproximation(timeStep));
        }
        if (name == getActionName<reactions::DetailedBalance>()) {
            return std::unique_ptr<TimeStepDependentAction>(detailedBalance(timeStep));
        }
        throw std::invalid_argument("Requested reaction scheduler " + name + " is not available.");
    }

    virtual std::unique_ptr<AddParticles> addParticles(const std::vector<Particle> &particles) const = 0;

    std::unique_ptr<AddParticles> addParticles(const Particle &particle) const {
        return addParticles(std::vector<Particle>{particle});
    }

    virtual std::unique_ptr<EulerBDIntegrator> eulerBDIntegrator(scalar timeStep) const = 0;

    virtual std::unique_ptr<MdgfrdIntegrator> mdgfrdIntegrator(scalar timeStep) const = 0;

    virtual std::unique_ptr<readdy::model::actions::CalculateForces> calculateForces() const = 0;

    virtual std::unique_ptr<CreateNeighborList> createNeighborList(scalar interactionDistance) const = 0;

    virtual std::unique_ptr<UpdateNeighborList> updateNeighborList() const = 0;

    virtual std::unique_ptr<ClearNeighborList> clearNeighborList() const = 0;

    virtual std::unique_ptr<EvaluateCompartments> evaluateCompartments() const = 0;

    virtual std::unique_ptr<reactions::UncontrolledApproximation>
    uncontrolledApproximation(scalar timeStep) const = 0;

    virtual std::unique_ptr<reactions::Gillespie>
    gillespie(scalar timeStep) const = 0;

    virtual std::unique_ptr<reactions::DetailedBalance>
    detailedBalance(scalar timeStep) const = 0;

    virtual std::unique_ptr<top::EvaluateTopologyReactions> evaluateTopologyReactions(scalar timeStep) const = 0;

};


NAMESPACE_END(actions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
