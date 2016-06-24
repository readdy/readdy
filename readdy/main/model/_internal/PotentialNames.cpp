/**
 * << detailed description >>
 *
 * @file PotentialNames.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */

#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/potentials/PotentialsOrder2.h>

namespace readdy {
    namespace model {
        namespace potentials {
            namespace _internal {
                // order 1
                const std::string PotentialName<CubePotential>::value = "Cube";

                // order 2
                const std::string PotentialName<HarmonicRepulsion>::value = "HarmonicRepulsion";
            }
        }
    }
}