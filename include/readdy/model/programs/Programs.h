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

#ifndef READDY_MAIN_PROGRAMS_H_H
#define READDY_MAIN_PROGRAMS_H_H

#include <readdy/model/programs/Program.h>
#include <readdy/model/Particle.h>
#include <type_traits>

#if READDY_OSX
#include <functional>
#endif

namespace readdy {
namespace model {
namespace programs {

class Test : public Program {

public:
    Test() : Program() {}
};

class AddParticle : public Program {

public:
    AddParticle() : Program() {}
};

class EulerBDIntegrator : public Program {
public:
    EulerBDIntegrator() : Program() {}
};

class CalculateForces : public Program {
public:
    CalculateForces() : Program() {}
};

class UpdateNeighborList : public Program {
public:
    enum Action {
        create, clear
    };

    UpdateNeighborList() : Program() {}

    void setAction(Action action) { UpdateNeighborList::action = action; }

protected:
    Action action = Action::create;
};

namespace reactions {
class UncontrolledApproximation : public Program {

public:
    using reaction_11 = std::function<model::Particle(const model::Particle &)>;
    using reaction_12 = std::function<void(const model::Particle &, model::Particle &, model::Particle &)>;
    using reaction_21 = std::function<model::Particle(const model::Particle &, const model::Particle &)>;
    using reaction_22 = std::function<void(const model::Particle &, const model::Particle &, model::Particle &,
                                           model::Particle &)>;

    UncontrolledApproximation() : Program() {}

    virtual void registerReactionScheme_11(const std::string &reactionName, reaction_11 fun) = 0;

    virtual void registerReactionScheme_12(const std::string &reactionName, reaction_12 fun) = 0;

    virtual void registerReactionScheme_21(const std::string &reactionName, reaction_21 fun) = 0;

    virtual void registerReactionScheme_22(const std::string &reactionName, reaction_22 fun) = 0;
};

class Gillespie : public Program {
public:
    Gillespie() : Program() {}
};

class GillespieParallel : public Program {
public:
    GillespieParallel() : Program() {}
};
struct NextSubvolumes : Program {
    NextSubvolumes() : Program() {}
};
}

class Compartments : public Program {
public:
    Compartments() : Program() {}

    virtual void registerCompartment(const std::function<bool(const readdy::model::Vec3)> characteristicFun) = 0;

    virtual void registerConversion(size_t compartmentIdx, std::string from, std::string to) = 0;
};

template<typename T>
const std::string getProgramName(typename std::enable_if<std::is_base_of<Test, T>::value>::type * = 0) {
    return "Test";
};

template<typename T>
const std::string getProgramName(typename std::enable_if<std::is_base_of<AddParticle, T>::value>::type * = 0) {
    return "AddParticle";
};

template<typename T>
const std::string getProgramName(typename std::enable_if<std::is_base_of<EulerBDIntegrator, T>::value>::type * = 0) {
    return "EulerBDIntegrator";
};

template<typename T>
const std::string getProgramName(typename std::enable_if<std::is_base_of<CalculateForces, T>::value>::type * = 0) {
    return "Calculate forces";
};

template<typename T>
const std::string getProgramName(typename std::enable_if<std::is_base_of<UpdateNeighborList, T>::value>::type * = 0) {
    return "Update neighbor list";
};

template<typename T>
const std::string
getProgramName(typename std::enable_if<std::is_base_of<reactions::UncontrolledApproximation, T>::value>::type * = 0) {
    return "UncontrolledApproximation";
};

template<typename T>
const std::string getProgramName(typename std::enable_if<std::is_base_of<reactions::Gillespie, T>::value>::type * = 0) {
    return "Gillespie";
};

template<typename T>
const std::string
getProgramName(typename std::enable_if<std::is_base_of<reactions::GillespieParallel, T>::value>::type * = 0) {
    return "GillespieParallel";
};

template<typename T>
const std::string
getProgramName(typename std::enable_if<std::is_base_of<reactions::NextSubvolumes, T>::value>::type * = 0) {
    return "NextSubvolumes";
};


template<typename T>
const std::string getProgramName(typename std::enable_if<std::is_base_of<Compartments, T>::value>::type * = 0) {
    return "Compartments";
};

}
}
}

#endif //READDY_MAIN_PROGRAMS_H_H
