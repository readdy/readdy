/**
 * Declaration of the base class of all order 1 potentials.
 * Subclasses have to implement calculateEnergy, calculateForce and calculateForceAndEnergy.
 * The first three methods take a modifiable reference and a particle's position. The last method is for replication
 * of the potential, so that it can be assigned to multiple particle types.
 *
 * @file PotentialOrder1.h
 * @brief Declaration of the base class of all order 1 potentials
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
    PotentialOrder1(const std::string &name) : Potential(name, 1) {}

    virtual double calculateEnergy(const Vec3 &position) const = 0;

    virtual void calculateForce(Vec3 &force, const Vec3 &position) const = 0;

    virtual void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const = 0;

    virtual void configureForType(const unsigned int type) {}

    virtual double getRelevantLengthScale() const noexcept = 0;

};
}
}
}
#endif //READDY_MAIN_POTENTIALORDER1_H
