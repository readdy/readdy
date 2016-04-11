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
            TestProgram() : Program("TestProgram") { }
            virtual void execute() = 0;
        };
    }
}

#endif //READDY2_MAIN_PROGRAMS_H_H
