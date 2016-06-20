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
        namespace programs {
            class Program {
            public:
                Program(const std::string &name) : name(name) { };

                virtual void execute() = 0;

                virtual ~Program() = default;

            protected:
                const std::string name;
            };
            namespace _internal {
                template<typename T>
                struct ProgramName {
                    static const std::string value;
                };
            }

            template<typename ProgramType>
            const std::string &getProgramName() {
                return readdy::model::programs::_internal::ProgramName<ProgramType>::value;
            }
        }
    }
}


#endif //READDY_MAIN_PROGRAM_H
