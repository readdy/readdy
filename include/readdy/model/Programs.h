//
// Created by clonker on 11.04.16.
//

#ifndef READDY2_MAIN_PROGRAMS_H_H
#define READDY2_MAIN_PROGRAMS_H_H

#include <readdy/model/Program.h>

namespace readdy {
    namespace model {
        class TestProgram : public Program {

        public:
            TestProgram() : Program() { }
            virtual ~TestProgram() = default;

            static std::string getName() {
                return "TestProgram";
            };
        };

        class AddParticleProgram : public Program {

        public:
            AddParticleProgram() : Program() {}
            virtual ~AddParticleProgram() = default;

            static std::string getName() {
                return "AddParticle";
            }
        };

        class DiffuseProgram : public Program {
        public:
            DiffuseProgram() : Program() {}
            virtual ~DiffuseProgram() = default;

            static std::string getName() {
                return "Diffuse";
            }
        };
    }
}

#endif //READDY2_MAIN_PROGRAMS_H_H
