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

    template<typename R, typename... Args>
    std::unique_ptr<R> createAction(Args &&... args) const {
        return std::unique_ptr<R>(get_dispatcher<R, Args...>::impl(this, std::forward<Args>(args)...));
    }

    std::vector<std::string> getAvailableActions() const {
        return {
                getActionName<AddParticles>(), getActionName<EulerBDIntegrator>(), getActionName<CalculateForces>(),
                getActionName<UpdateNeighborList>(), getActionName<reactions::UncontrolledApproximation>(),
                getActionName<reactions::Gillespie>(), getActionName<reactions::GillespieParallel>(),
                getActionName<reactions::NextSubvolumes>()
        };
    }

    /*
     * Convenience stuff
     */
    std::unique_ptr<TimeStepDependentAction> createIntegrator(const std::string& name, double timeStep) {
        if(name == getActionName<EulerBDIntegrator>()) {
            return std::unique_ptr<TimeStepDependentAction>(createEulerBDIntegrator(timeStep));
        }
        log::critical("Requested integrator \"{}\" is not available, returning nullptr", name);
        return nullptr;
    }

    std::unique_ptr<TimeStepDependentAction> createReactionScheduler(const std::string& name, double timeStep) {
        if(name == getActionName<reactions::Gillespie>()) {
            return std::unique_ptr<TimeStepDependentAction>(createGillespie(timeStep));
        } else if(name == getActionName<reactions::GillespieParallel>()) {
            return std::unique_ptr<TimeStepDependentAction>(createGillespieParallel(timeStep));
        } else if(name == getActionName<reactions::NextSubvolumes>()) {
            return std::unique_ptr<TimeStepDependentAction>(createNextSubvolumes(timeStep));
        } else if(name == getActionName<reactions::UncontrolledApproximation>()) {
            return std::unique_ptr<TimeStepDependentAction>(createUncontrolledApproximation(timeStep));
        }
        log::critical("Requested reaction scheduler \"{}\" is not available, returning nullptr", name);
        return nullptr;
    }

protected:
    
    virtual AddParticles *createAddParticles(const std::vector<Particle> &particles) const = 0;
    AddParticles *createAddParticles(const Particle &particle) const {
        return createAddParticles(std::vector<Particle>{particle});
    };

    virtual EulerBDIntegrator *createEulerBDIntegrator(double timeStep) const = 0;

    virtual CalculateForces *createCalculateForces() const = 0;

    virtual UpdateNeighborList *createUpdateNeighborList(UpdateNeighborList::Operation, double skinSize) const = 0;
    UpdateNeighborList *createUpdateNeighborList(UpdateNeighborList::Operation op) const {
        return createUpdateNeighborList(op, -1);
    };
    UpdateNeighborList *createUpdateNeighborList() const {
        return createUpdateNeighborList(UpdateNeighborList::Operation::create, -1);
    };

    virtual EvaluateCompartments *createEvaluateCompartments() const = 0;

    virtual reactions::UncontrolledApproximation *createUncontrolledApproximation(double timeStep) const = 0;

    virtual reactions::Gillespie *createGillespie(double timeStep) const = 0;

    virtual reactions::GillespieParallel *createGillespieParallel(double timeStep) const = 0;

    virtual reactions::NextSubvolumes *createNextSubvolumes(double timeStep) const = 0;

    template<typename T, typename... Args>
    struct get_dispatcher;

    template<typename T, typename... Args>
    struct get_dispatcher {
        static T *impl(const ActionFactory *self, Args &&... args) {
            // this only invokes the normal constructor
            return new T(std::forward<Args>(args)...);
        };
    };
};

READDY_CREATE_FACTORY_DISPATCHER(ActionFactory, AddParticles)

READDY_CREATE_FACTORY_DISPATCHER(ActionFactory, EulerBDIntegrator)

READDY_CREATE_FACTORY_DISPATCHER(ActionFactory, CalculateForces)

READDY_CREATE_FACTORY_DISPATCHER(ActionFactory, UpdateNeighborList)

READDY_CREATE_FACTORY_DISPATCHER(ActionFactory, EvaluateCompartments)

READDY_CREATE_FACTORY_DISPATCHER2(ActionFactory, reactions, UncontrolledApproximation)

READDY_CREATE_FACTORY_DISPATCHER2(ActionFactory, reactions, Gillespie)

READDY_CREATE_FACTORY_DISPATCHER2(ActionFactory, reactions, GillespieParallel)

READDY_CREATE_FACTORY_DISPATCHER2(ActionFactory, reactions, NextSubvolumes)

NAMESPACE_END(actions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
