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
 * This file contains the declaration of the program factory. Internally, the factory is simply a map of
 * string -> std::function<Program*()>, which then can get called.
 *
 * @file ProgramFactory.h
 * @brief Declaration of the program factory.
 * @author clonker
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
                getActionName<UpdateNeighborList>(), getActionName<reactions::UncontrolledApproximation>(),
                getActionName<reactions::Gillespie>(), getActionName<reactions::DetailedBalance>(),
                /*getActionName<reactions::GillespieParallel>(),*/
                /*getActionName<reactions::NextSubvolumes>(),*/ getActionName<top::EvaluateTopologyReactions>()
        };
    }

    /*
     * Convenience stuff
     */
    std::unique_ptr<TimeStepDependentAction> createIntegrator(const std::string& name, scalar timeStep) {
        if(name == getActionName<EulerBDIntegrator>()) {
            return std::unique_ptr<TimeStepDependentAction>(eulerBDIntegrator(timeStep));
        }
        log::critical("Requested integrator \"{}\" is not available, returning nullptr", name);
        return nullptr;
    }

    std::unique_ptr<TimeStepDependentAction> createReactionScheduler(const std::string& name, scalar timeStep) {
        if(name == getActionName<reactions::Gillespie>()) {
            return std::unique_ptr<TimeStepDependentAction>(gillespie(timeStep));
        }
        if(name == getActionName<reactions::UncontrolledApproximation>()) {
            return std::unique_ptr<TimeStepDependentAction>(uncontrolledApproximation(timeStep));
        }
        if(name==getActionName<reactions::DetailedBalance>()) {
            return std::unique_ptr<TimeStepDependentAction>(detailedBalance(timeStep));
        }
        log::critical("Requested reaction scheduler \"{}\" is not available, returning nullptr", name);
        return nullptr;
    }

    virtual std::unique_ptr<AddParticles> addParticles(const std::vector<Particle> &particles) const = 0;
    std::unique_ptr<AddParticles> addParticles(const Particle &particle) const {
        return addParticles(std::vector<Particle>{particle});
    }

    virtual std::unique_ptr<EulerBDIntegrator> eulerBDIntegrator(scalar timeStep) const = 0;

    virtual std::unique_ptr<CalculateForces> calculateForces() const = 0;

    virtual std::unique_ptr<UpdateNeighborList> updateNeighborList(UpdateNeighborList::Operation operation,
                                                                   scalar skinSize) const = 0;

    std::unique_ptr<UpdateNeighborList> updateNeighborList(UpdateNeighborList::Operation op) const {
        return updateNeighborList(op, 0);
    };
    std::unique_ptr<UpdateNeighborList> updateNeighborList() const {
        return updateNeighborList(UpdateNeighborList::Operation::init, 0);
    };

    virtual std::unique_ptr<EvaluateCompartments> evaluateCompartments() const = 0;

    virtual std::unique_ptr<reactions::UncontrolledApproximation> uncontrolledApproximation(scalar timeStep) const = 0;

    virtual std::unique_ptr<reactions::Gillespie> gillespie(scalar timeStep) const = 0;

    virtual std::unique_ptr<reactions::DetailedBalance> detailedBalance(scalar timeStep) const = 0;

    virtual std::unique_ptr<top::EvaluateTopologyReactions> evaluateTopologyReactions(scalar timeStep) const = 0;

};


NAMESPACE_END(actions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
