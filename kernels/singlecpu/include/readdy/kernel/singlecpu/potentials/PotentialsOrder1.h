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

namespace readdy {
    namespace kernel {
        namespace singlecpu {

            class SingleCPUKernel;
            namespace potentials {

                class SingleCPUCubePotential : public readdy::model::potentials::CubePotential {
                public:
                    SingleCPUCubePotential(const SingleCPUKernel *const kernel);

                    virtual double calculateEnergy(const readdy::model::Vec3 &position) override;

                    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) override;

                    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &position) override;

                    virtual potentials::SingleCPUCubePotential *replicate() const override;


                };

            }
        }
    }

}

#endif //READDY_MAIN_P1CUBE_H
