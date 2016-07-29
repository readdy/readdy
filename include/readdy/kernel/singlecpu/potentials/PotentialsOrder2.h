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
#include <readdy/model/Kernel.h>

namespace readdy {

    namespace kernel {
        namespace singlecpu {

            namespace potentials {
                class SingleCPUHarmonicRepulsion : public readdy::model::potentials::HarmonicRepulsion{
                    using vec_t = readdy::model::Vec3;

                public:
                    SingleCPUHarmonicRepulsion(const readdy::model::Kernel *const kernel);

                    virtual double calculateEnergy(const vec_t& x_ij) override;
                    virtual void calculateForce(vec_t &force, const vec_t &x_ij) override;
                    virtual void calculateForceAndEnergy(vec_t &force, double &energy, const vec_t &x_ij) override;

                    virtual double getCutoffRadius() const override;


                    virtual potentials::SingleCPUHarmonicRepulsion *replicate() const override;


                };

                class SingleCPUWeakInteractionPiecewiseHarmonic : public readdy::model::potentials::WeakInteractionPiecewiseHarmonic {
                    using vec_t = readdy::model::Vec3;
                public:
                    SingleCPUWeakInteractionPiecewiseHarmonic(const readdy::model::Kernel *const kernel);

                    virtual double calculateEnergy(const readdy::model::Vec3 &x_ij) override;
                    virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &x_ij) override;
                    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &x_ij) override;

                    virtual potentials::SingleCPUWeakInteractionPiecewiseHarmonic *replicate() const override;

                    virtual double getCutoffRadius() const override;


                };

            }
        }
    }
}
#endif //READDY_MAIN_POTENTIALSORDER2_H_H
