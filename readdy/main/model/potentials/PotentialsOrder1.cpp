/**
 * << detailed description >>
 *
 * @file PotentialsOrder1.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 31.05.16
 */

#include <readdy/model/potentials/PotentialsOrder1.h>

namespace readdy {
    namespace model {
        namespace potentials {
            P1Cube::P1Cube(const unsigned int id)  : PotentialOrder1(id, getPotentialName<P1Cube>()) { }

            namespace _internal {
                const std::string PotentialName<P1Cube>::value = "P1Cube";
            }
        }
    }
}

