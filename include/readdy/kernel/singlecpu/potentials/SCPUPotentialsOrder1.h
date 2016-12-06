/**
 * << detailed description >>
 *
 * @file P1Cube.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#ifndef READDY_MAIN_P1CUBE_H
#define READDY_MAIN_P1CUBE_H

#include <string>
#include <readdy/model/potentials/PotentialOrder1.h>
#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/Kernel.h>

namespace readdy {
namespace kernel {
namespace scpu {

namespace potentials {

class SCPUCubePotential : public readdy::model::potentials::CubePotential {
public:
    SCPUCubePotential(const readdy::model::Kernel *const kernel);

    virtual double calculateEnergy(const readdy::model::Vec3 &position) const override;

    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) const override;

    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy,
                                         const readdy::model::Vec3 &position) const override;

};

class SCPUSpherePotential : public readdy::model::potentials::SpherePotential {
public:
    SCPUSpherePotential(const readdy::model::Kernel * const kernel);

    virtual double calculateEnergy(const readdy::model::Vec3 &position) const override;

    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) const override;

    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &position) const override;

};

}
}
}

}

#endif //READDY_MAIN_P1CUBE_H
