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

                public:
                    HarmonicRepulsion(const SingleCPUKernel * const);

                    virtual double calculateEnergy(const size_t &i, const size_t &j) override;
                    virtual void calculateForce(readdy::model::Vec3 &force, const size_t &i, const size_t &j) override;
                    virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const size_t &i, const size_t &j) override;


                };

            }
        }
    }

    namespace model {
        namespace potentials {
            namespace _internal {
                namespace pot = readdy::kernel::singlecpu::potentials;
                template<> struct PotentialName<pot::HarmonicRepulsion> { static const std::string value; };
            }
        }
    }
}
#endif //READDY_MAIN_POTENTIALSORDER2_H_H
