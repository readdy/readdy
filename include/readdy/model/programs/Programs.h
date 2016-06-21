//
// Created by clonker on 11.04.16.
//

#ifndef READDY_MAIN_PROGRAMS_H_H
#define READDY_MAIN_PROGRAMS_H_H

#include <readdy/model/programs/Program.h>
#include <readdy/common/Types.h>
#include <readdy/model/Particle.h>

namespace readdy {
    namespace model {
        namespace programs {
            class TestProgram : public Program {

            public:
                TestProgram() : Program(getProgramName<TestProgram>()) { }
            };

            class AddParticleProgram : public Program {

            public:
                AddParticleProgram() : Program(getProgramName<AddParticleProgram>()) { }
            };

            class DiffuseProgram : public Program {
            public:
                DiffuseProgram() : Program(getProgramName<DiffuseProgram>()) { }
            };

            class UpdateStateModelProgram : public Program {
            public:
                UpdateStateModelProgram() : Program(getProgramName<UpdateStateModelProgram>()) { }
                void configure(const readdy::model::time_step_type& t) {
                    curr_t = t;
                }
            protected:
                readdy::model::time_step_type curr_t;
            };

            class DefaultReactionProgram : public Program {

            public:
                using reaction_11 = std::function<model::Particle(const model::Particle&)>;
                using reaction_12 = std::function<void(const model::Particle&, model::Particle&, model::Particle&)>;
                using reaction_21 = std::function<model::Particle(const model::Particle&, const model::Particle&)>;
                using reaction_22 = std::function<void(const model::Particle&, const model::Particle&, model::Particle&, model::Particle&)>;

                DefaultReactionProgram() : Program(getProgramName<DefaultReactionProgram>()) { }
                virtual void configure() = 0;
                virtual void registerReactionScheme_11(const std::string& reactionName, reaction_11 fun) = 0;
                virtual void registerReactionScheme_12(const std::string& reactionName, reaction_12 fun) = 0;
                virtual void registerReactionScheme_21(const std::string& reactionName, reaction_21 fun) = 0;
                virtual void registerReactionScheme_22(const std::string& reactionName, reaction_22 fun) = 0;
            };

            namespace _internal {
                template<> struct ProgramName<TestProgram> { static const std::string value; };
                template<> struct ProgramName<AddParticleProgram> { static const std::string value; };
                template<> struct ProgramName<DiffuseProgram> { static const std::string value; };
                template<> struct ProgramName<UpdateStateModelProgram> { static const std::string value; };
                template<> struct ProgramName<DefaultReactionProgram> { static const std::string value; };
            }
        }
    }
}

#endif //READDY_MAIN_PROGRAMS_H_H
