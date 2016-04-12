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
            Program(const std::string name);
            virtual ~Program();
            Program(Program&& rhs);
            Program& operator=(Program&& rhs);

            const std::string getName() const;
            virtual void execute() = 0;

        private:
            struct Impl;
            std::unique_ptr<Impl> impl_ptr;
        };
    }
}


#endif //READDY2_MAIN_PROGRAM_H
