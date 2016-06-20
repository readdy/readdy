//
// Created by clonker on 11.04.16.
//

#ifndef READDY_MAIN_PROGRAMS_H_H
#define READDY_MAIN_PROGRAMS_H_H

#include <readdy/model/programs/Program.h>
#include <readdy/common/Types.h>

namespace readdy {
    namespace model {
        namespace programs {
            class TestProgram : public Program {

            public:
                TestProgram() : Program(getProgramName<TestProgram>()) { }

                virtual ~TestProgram() = default;
            };

            class AddParticleProgram : public Program {

            public:
                AddParticleProgram() : Program(getProgramName<AddParticleProgram>()) { }

                virtual ~AddParticleProgram() = default;
            };

            class DiffuseProgram : public Program {
            public:
                DiffuseProgram() : Program(getProgramName<DiffuseProgram>()) { }

                virtual ~DiffuseProgram() = default;
            };

            class UpdateStateModelProgram : public Program {
            public:
                UpdateStateModelProgram() : Program(getProgramName<UpdateStateModelProgram>()) { }

                virtual ~UpdateStateModelProgram() = default;

                void execute_t(readdy::model::time_step_type t) {
                        setCurrentTimeStep(t);
                        execute();
                }

                virtual void setCurrentTimeStep(readdy::model::time_step_type t) = 0;
            };

            namespace _internal {
                template<>
                struct ProgramName<TestProgram> {
                    static const std::string value;
                };
                template<>
                struct ProgramName<AddParticleProgram> {
                    static const std::string value;
                };
                template<>
                struct ProgramName<DiffuseProgram> {
                    static const std::string value;
                };
                template<>
                struct ProgramName<UpdateStateModelProgram> {
                    static const std::string value;
                };
            }
        }
    }
}

#endif //READDY_MAIN_PROGRAMS_H_H
