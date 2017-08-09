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
 * This files contains a selection of possible programs that can be executed on a kernel:
 *   - TestProgram: Program that has no other operation than printing something to the log.
 *   - AddParticleProgram: A program with which particles can be added.
 *   - EulerBDIntegrator: A program that propagates the particles through the system. The update model program should be
 *                     called beforehand, such that forces are available.
 *   - UpdateNeighborList: A program that creates neighbor lists.
 *   - CalculateForces: A program that calculates forces for later use in, e.g., integration schemes.
 *   - DefaultReactionProgram: A program that executes the default reaction scheme.
 *   - CompartmentConversion: Perform instantaneous particle conversions depending on the particles' position.
 *
 * Further, specializations of ProgramName<T> are declared.
 *
 * @file Programs.h
 * @brief Declaration of all globally available programs.
 * @author clonker
 * @date 11.04.16
 * @todo provide more detailed descriptions for some of the programs
 */
#pragma once

#include <type_traits>
#include <readdy/model/Particle.h>
#include <readdy/model/actions/Action.h>
#include <readdy/model/observables/Observable.h>

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

    void perform() override;

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
        create, clear
    };

    explicit UpdateNeighborList(Operation operation = Operation::create, scalar skinSize = -1);

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

class GillespieParallel : public TimeStepDependentAction {
public:
    explicit GillespieParallel(scalar timeStep);
};

struct NextSubvolumes : public TimeStepDependentAction {
    explicit NextSubvolumes(scalar timeStep);
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
getActionName(typename std::enable_if<std::is_base_of<reactions::GillespieParallel, T>::value>::type * = 0) {
    return "GillespieParallel";
}

template<typename T>
const std::string
getActionName(typename std::enable_if<std::is_base_of<reactions::NextSubvolumes, T>::value>::type * = 0) {
    return "NextSubvolumes";
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
