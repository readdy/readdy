/**
 * This files contains a selection of possible programs that can be executed on a kernel:
 *   - TestProgram: Program that has no other operation than printing something to the log.
 *   - AddParticleProgram: A program with which particles can be added.
 *   - EulerBDIntegrator: A program that propagates the particles through the system. The update model program should be
 *                     called beforehand, such that forces are available.
 *   - UpdateNeighborList: A program that creates neighbor lists.
 *   - CalculateForces: A program that calculates forces for later use in, e.g., integration schemes.
 *   - DefaultReactionProgram: A program that executes the default reaction scheme.
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
#include <readdy/common/Types.h>
#include <readdy/model/Particle.h>

namespace readdy {
    namespace model {
        namespace programs {
            class Test : public Program {

            public:
                Test() : Program(getProgramName<Test>()) { }
            };

            class AddParticle : public Program {

            public:
                AddParticle() : Program(getProgramName<AddParticle>()) { }
            };

            class EulerBDIntegrator : public Program {
            public:
                EulerBDIntegrator() : Program(getProgramName<EulerBDIntegrator>()) { }
            };

            class CalculateForces : public Program {
            public:
                CalculateForces() : Program(getProgramName<CalculateForces>()) {}
            };

            class UpdateNeighborList : public Program {
            public:
                UpdateNeighborList() : Program(getProgramName<UpdateNeighborList>()) {}
            };

            namespace reactions {
                class UncontrolledApproximation : public Program {

                public:
                    using reaction_11 = std::function<model::Particle(const model::Particle&)>;
                    using reaction_12 = std::function<void(const model::Particle&, model::Particle&, model::Particle&)>;
                    using reaction_21 = std::function<model::Particle(const model::Particle&, const model::Particle&)>;
                    using reaction_22 = std::function<void(const model::Particle&, const model::Particle&, model::Particle&, model::Particle&)>;

                    UncontrolledApproximation() : Program(getProgramName<UncontrolledApproximation>()) { }
                    virtual void registerReactionScheme_11(const std::string& reactionName, reaction_11 fun) = 0;
                    virtual void registerReactionScheme_12(const std::string& reactionName, reaction_12 fun) = 0;
                    virtual void registerReactionScheme_21(const std::string& reactionName, reaction_21 fun) = 0;
                    virtual void registerReactionScheme_22(const std::string& reactionName, reaction_22 fun) = 0;
                };



            }

            namespace _internal {
                template<> struct ProgramName<Test> { static const std::string value; };
                template<> struct ProgramName<AddParticle> { static const std::string value; };
                template<> struct ProgramName<EulerBDIntegrator> { static const std::string value; };
                template<> struct ProgramName<CalculateForces> { static const std::string value; };
                template<> struct ProgramName<UpdateNeighborList> { static const std::string value; };

                template<> struct ProgramName<reactions::UncontrolledApproximation> { static const std::string value; };
            }
        }
    }
}

#endif //READDY_MAIN_PROGRAMS_H_H
