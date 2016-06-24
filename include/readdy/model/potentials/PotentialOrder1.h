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
                virtual double calculateEnergy(const Vec3& position) = 0;
                virtual void calculateForce(Vec3 &force, const Vec3& position) = 0;
                virtual void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3& position) = 0;

                virtual void configureForType(const unsigned int &type) {}

                virtual PotentialOrder1 *replicate() const override = 0;


            };
        }
    }
}
#endif //READDY_MAIN_POTENTIALORDER1_H
