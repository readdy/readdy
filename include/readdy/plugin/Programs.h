//
// Created by clonker on 11.04.16.
//

#ifndef READDY2_MAIN_PROGRAMS_H_H
#define READDY2_MAIN_PROGRAMS_H_H

#include <readdy/plugin/Program.h>

namespace readdy {
    namespace plugin {
        class TestProgram : public Program {

        public:
            TestProgram() : Program() { }

            static const std::string getName() {
                return "TestProgram";
            };
        };
    }
}

#endif //READDY2_MAIN_PROGRAMS_H_H
