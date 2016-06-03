/**
 * << detailed description >>
 *
 * @file PotentialFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 31.05.16
 */

#ifndef READDY_MAIN_POTENTIALFACTORY_H
#define READDY_MAIN_POTENTIALFACTORY_H

#include <unordered_map>
#include <functional>

#include <readdy/model/potentials/Potential.h>

namespace readdy {
    namespace model {

        class Kernel;
        namespace potentials {
            class PotentialFactory {
            public:
                PotentialFactory(Kernel *const kernel);

            protected:
                std::unordered_map<std::string, std::function<Potential *()>> factory;
                Kernel *const kernel;
            };

        }
    }
}

#endif //READDY_MAIN_POTENTIALFACTORY_H
