//
// Created by clonker on 08.04.16.
//

#ifndef READDY2_MAIN_PROGRAMFACTORY_H
#define READDY2_MAIN_PROGRAMFACTORY_H

#include <type_traits>
#include <readdy/plugin/Program.h>

namespace readdy {
    namespace plugin {
        template<class T>
        class ProgramFactory {
            /*static_assert(std::is_base_of<ProgramFactory, T>::value,
                          "The given template parameter did not refer to a class that has ProgramFactory as its base.");*/

            static Program *createProgram(const std::string name) {
                return T::createProgram(name);
            }
        };
    }
}

#endif //READDY2_MAIN_PROGRAMFACTORY_H
