//
// Created by clonker on 08.04.16.
//

#ifndef READDY2_MAIN_PROGRAMFACTORY_H
#define READDY2_MAIN_PROGRAMFACTORY_H

#include <type_traits>
#include <readdy/plugin/Program.h>

namespace readdy {
    namespace plugin {
        template<class FactorySpecialization>
        class ProgramFactory {
            static std::shared_ptr<Program> createProgram(const std::string name) {
                return FactorySpecialization::createProgram(name);
            }
        };
    }
}

#endif //READDY2_MAIN_PROGRAMFACTORY_H
