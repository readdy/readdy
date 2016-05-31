//
// Created by clonker on 08.04.16.
//

#ifndef READDY_MAIN_PROGRAM_H
#define READDY_MAIN_PROGRAM_H

#include <memory>

#include <boost/predef.h>
#if BOOST_OS_MACOS
#include <string>
#endif

namespace readdy {
    namespace model {
        class Program {
        public:
            virtual void execute() = 0;
            virtual ~Program() = default;
        };
    }
}


#endif //READDY_MAIN_PROGRAM_H
