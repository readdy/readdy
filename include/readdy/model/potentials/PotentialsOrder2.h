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
        class Kernel;
        namespace potentials {

            class HarmonicRepulsion : public PotentialOrder2 {

            public:
                HarmonicRepulsion(const Kernel * const kernel);

                virtual HarmonicRepulsion *replicate() const override = 0;

                double getSumOfParticleRadii() const;
                double getSumOfParticleRadiiSquared() const;
                double getForceConstant() const;
                void setForceConstant(double forceConstant);

                virtual void configureForTypes(unsigned int type1, unsigned int type2) override;

            protected:
                const Kernel * const kernel;
                double sumOfParticleRadii;
                double sumOfParticleRadiiSquared;
                double forceConstant;

            };

            namespace _internal {
                template<> struct PotentialName<HarmonicRepulsion> { static const std::string value; };
            }
        }
    }
}

#endif //READDY_MAIN_POTENTIALSORDER2_H
