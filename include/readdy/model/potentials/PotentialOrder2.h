/**
 * << detailed description >>
 *
 * @file PotentialOrder2.h
 * @brief << brief description >>
 * @author clonker
 * @date 31.05.16
 */

#ifndef READDY_MAIN_POTENTIALORDER2_H
#define READDY_MAIN_POTENTIALORDER2_H

#include "Potential.h"

namespace readdy {
    namespace model {
        namespace potentials {
            class PotentialOrder2 : public Potential {

            public:
                PotentialOrder2(const unsigned int id, const std::string &name) : Potential(id, name, 2) { }
            };
        }
    }
}
#endif //READDY_MAIN_POTENTIALORDER2_H
