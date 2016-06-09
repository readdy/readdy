/**
 * << detailed description >>
 *
 * @file PotentialOrder1.h
 * @brief << brief description >>
 * @author clonker
 * @date 31.05.16
 */

#ifndef READDY_MAIN_POTENTIALORDER1_H
#define READDY_MAIN_POTENTIALORDER1_H

#include <readdy/model/potentials/Potential.h>

namespace readdy {
    namespace model {
        namespace potentials {
            class PotentialOrder1 : public Potential {

            public:
                PotentialOrder1(const std::string &name) : Potential(name, 1) { }
                virtual double calculateEnergy(size_t& i) = 0;
                virtual void calculateForce(Vec3 &force, size_t& i) = 0;
                virtual void calculateForceAndEnergy(Vec3 &force, double &energy, size_t& i) = 0;
            };
        }
    }
}
#endif //READDY_MAIN_POTENTIALORDER1_H
