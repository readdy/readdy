/**
 * << detailed description >>
 *
 * @file PotentialsOrder2.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#ifndef READDY_MAIN_POTENTIALSORDER2_H
#define READDY_MAIN_POTENTIALSORDER2_H

#include <readdy/model/KernelContext.h>
#include "PotentialOrder2.h"

namespace readdy {
    namespace model {
        namespace potentials {

            template<typename KernelType>
            class HarmonicRepulsion : public PotentialOrder2 {

            public:
                HarmonicRepulsion(const KernelType * const kernel) :
                        PotentialOrder2(_internal::PotentialName<HarmonicRepulsion<KernelType>>::value),
                        kernel(kernel) { }

                double getSumOfParticleRadii() const {
                    return sumOfParticleRadii;
                }


                double getSumOfParticleRadiiSquared() const {
                    return sumOfParticleRadiiSquared;
                }


                double getForceConstant() const {
                    return forceConstant;
                }

                void setForceConstant(double forceConstant) {
                    HarmonicRepulsion::forceConstant = forceConstant;
                }
                virtual void configureForTypes(unsigned int type1, unsigned int type2) override {
                    readdy::model::KernelContext ctx = kernel->getKernelContext();
                    auto r1 = ctx.getParticleRadius(type1);
                    auto r2 = ctx.getParticleRadius(type2);
                    sumOfParticleRadii = r1+r2;
                    sumOfParticleRadiiSquared = sumOfParticleRadii * sumOfParticleRadii;
                    // todo
                }

            protected:
                const KernelType * const kernel;
                double sumOfParticleRadii;
                double sumOfParticleRadiiSquared;
                double forceConstant;


            };

            namespace _internal {
                template<typename KernelType>
                struct PotentialName<HarmonicRepulsion<KernelType>> { static const std::string value; };
            }
        }
    }
}

#include "_impl/PotentialsOrder2Impl.h"

#endif //READDY_MAIN_POTENTIALSORDER2_H
