/**
 * << detailed description >>
 *
 * @file P1Cube.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#include "PotentialsOrder1.h"

namespace readdy {

    namespace kernel {
        namespace singlecpu {
            namespace potentials {
                P1Cube::P1Cube() : PotentialOrder1(getPotentialName<P1Cube>()) {

                }
            }
        }
    }


    namespace model {
        namespace potentials {
            namespace _internal {
                const std::string PotentialName<pot::P1Cube>::value = "P1Cube";
            }
        }
    }
}