/**
 * << detailed description >>
 *
 * @file Potentials.h
 * @brief << brief description >>
 * @author clonker
 * @date 31.05.16
 */

#ifndef READDY_MAIN_POTENTIALS_H
#define READDY_MAIN_POTENTIALS_H

#include "PotentialOrder1.h"

namespace readdy {
    namespace model {
        namespace potentials {

            template<typename PotentialType>
            const std::string& getPotentialName() {
                return _internal::PotentialName<PotentialType>::value;
            }

            class P1Cube : public PotentialOrder1 {
            public:
                P1Cube(const unsigned int id);
            };

            namespace _internal {
                template<> struct PotentialName<P1Cube> { static const std::string value; };
            }
        }
    }
}

#endif //READDY_MAIN_POTENTIALS_H
