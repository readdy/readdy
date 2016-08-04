/**
 * Declaration of the base class for all order 2 potentials. They basically have calculateForce, calculateEnergy and
 * calculateForceAndEnergy methods, which take a modifiable reference argument and the difference vector x_ij between
 * two particles.
 * Further, subclasses have to implement getCutoffRadius so that the neighbor list can be created more efficiently.
 *
 * @file PotentialOrder2.h
 * @brief Declaration of the base class for all order 2 potentials.
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
                virtual double calculateEnergy(const Vec3& x_ij) const = 0;
                virtual void calculateForce(Vec3 &force, const Vec3& x_ij) const = 0;
                virtual void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3& x_ij) const = 0;
                virtual void configureForTypes(unsigned int type1, unsigned int type2) {}

                virtual PotentialOrder2 *replicate() const override = 0;

                virtual double getCutoffRadius() const = 0;

            };
        }
    }
}
#endif //READDY_MAIN_POTENTIALORDER2_H
