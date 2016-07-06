/**
 * This header contains the declarations of order 2 potentials. Currently:
 *   - Harmonic repulsion
 *   - Weak interaction piecewise harmonic
 *
 * @file PotentialsOrder2.h
 * @brief Contains the declaration of order 2 potentials.
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
                double forceConstant = 0;

            };

            class WeakInteractionPiecewiseHarmonic : public PotentialOrder2 {

            public:
                WeakInteractionPiecewiseHarmonic(const Kernel* const kernel);
                virtual WeakInteractionPiecewiseHarmonic *replicate() const override = 0;

                virtual void configureForTypes(unsigned int type1, unsigned int type2) override;


                void setDesiredParticleDistance(double desiredParticleDistance);
                void setForceConstant(double forceConstant);
                void setDepthAtDesiredDistance(double depthAtDesiredDistance);
                void setNoInteractionDistance(double noInteractionDistance);

            protected:
                const Kernel* const kernel;
                double desiredParticleDistance;
                double forceConstant;
                double depthAtDesiredDistance;
                double noInteractionDistance;
            };

            namespace _internal {
                template<> struct PotentialName<HarmonicRepulsion> { static const std::string value; };
                template<> struct PotentialName<WeakInteractionPiecewiseHarmonic> { static const std::string value; };
            }
        }
    }
}

#endif //READDY_MAIN_POTENTIALSORDER2_H
