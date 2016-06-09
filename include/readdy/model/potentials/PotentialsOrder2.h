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

#include "PotentialOrder2.h"

namespace readdy {
    namespace model {
        namespace potentials {

            template<typename KernelType>
            class HarmonicRepulsion : public PotentialOrder2 {

            public:
                HarmonicRepulsion(const std::string& name, const KernelType * const kernel) : PotentialOrder2(name), kernel(kernel) { }

                double getSumOfParticleRadii() const {
                    return sumOfParticleRadii;
                }

                void setSumOfParticleRadii(double sumOfParticleRadii) {
                    HarmonicRepulsion::sumOfParticleRadii = sumOfParticleRadii;
                }

                double getSumOfParticleRadiiSquared() const {
                    return sumOfParticleRadiiSquared;
                }

                void setSumOfParticleRadiiSquared(double sumOfParticleRadiiSquared) {
                    HarmonicRepulsion::sumOfParticleRadiiSquared = sumOfParticleRadiiSquared;
                }

                double getForceConstant() const {
                    return forceConstant;
                }

                void setForceConstant(double forceConstant) {
                    HarmonicRepulsion::forceConstant = forceConstant;
                }

            protected:
                const KernelType * const kernel;
                double sumOfParticleRadii;
                double sumOfParticleRadiiSquared;
                double forceConstant;
            };

        }
    }
}
#endif //READDY_MAIN_POTENTIALSORDER2_H
