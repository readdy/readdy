//
// Created by clonker on 08.04.16.
//

#ifndef READDY2_MAIN_PROGRAM_H
#define READDY2_MAIN_PROGRAM_H

#include <memory>

#include <boost/predef.h>
#ifdef BOOST_OS_MACOS
#include <string>
#endif

namespace readdy {
    namespace plugin {
        class Program {
        public:
            virtual void execute() = 0;
            virtual ~Program() = default;
        };
    }
}


#endif //READDY2_MAIN_PROGRAM_H
