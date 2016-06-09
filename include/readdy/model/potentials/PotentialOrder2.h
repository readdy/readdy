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
                PotentialOrder2(const std::string &name) : Potential(name, 2) { }
                virtual double calculateEnergy(const size_t& i, const size_t& j) = 0;
                virtual void calculateForce(Vec3 &force, const size_t& i, const size_t& j) = 0;
                virtual void calculateForceAndEnergy(Vec3 &force, double &energy, const size_t& i, const size_t& j) = 0;

            };
        }
    }
}
#endif //READDY_MAIN_POTENTIALORDER2_H
