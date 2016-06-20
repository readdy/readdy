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
                class HarmonicRepulsion : public readdy::model::potentials::HarmonicRepulsion<SingleCPUKernel>{
                    using vec_t = readdy::model::Vec3;

                public:
                    HarmonicRepulsion(const SingleCPUKernel *const kernel);

                    virtual double calculateEnergy(const vec_t& x_i, const vec_t& x_j) override;
                    virtual void calculateForce(vec_t &force, const vec_t &x_i, const vec_t &x_j) override;
                    virtual void calculateForceAndEnergy(vec_t &force, double &energy, const vec_t &x_i, const vec_t &x_j) override;

                    virtual potentials::HarmonicRepulsion *replicate() const override;


                };

            }
        }
    }
}
#endif //READDY_MAIN_POTENTIALSORDER2_H_H
