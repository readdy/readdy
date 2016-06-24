/**
 * << detailed description >>
 *
 * @file PotentialsOrder2.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#ifndef READDY_MAIN_POTENTIALSORDER2_H_H
#define READDY_MAIN_POTENTIALSORDER2_H_H

#include <readdy/model/potentials/PotentialsOrder2.h>

namespace readdy {

    namespace kernel {
        namespace singlecpu {
            class SingleCPUKernel;

            namespace potentials {
                class SingleCPUHarmonicRepulsion : public readdy::model::potentials::HarmonicRepulsion{
                    using vec_t = readdy::model::Vec3;

                public:
                    SingleCPUHarmonicRepulsion(const SingleCPUKernel *const kernel);

                    virtual double calculateEnergy(const vec_t& x_ij) override;
                    virtual void calculateForce(vec_t &force, const vec_t &x_ij) override;
                    virtual void calculateForceAndEnergy(vec_t &force, double &energy, const vec_t &x_ij) override;

                    virtual potentials::SingleCPUHarmonicRepulsion *replicate() const override;


                };

            }
        }
    }
}
#endif //READDY_MAIN_POTENTIALSORDER2_H_H
