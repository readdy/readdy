/**
 * << detailed description >>
 *
 * @file PotentialsOrder1.h
 * @brief << brief description >>
 * @author clonker
 * @date 15.06.16
 */

#include "PotentialOrder1.h"

#ifndef READDY_MAIN_POTENTIALSORDER1_H
#define READDY_MAIN_POTENTIALSORDER1_H

namespace readdy {
    namespace model {
        class Kernel;
        namespace potentials {
            class CubePotential : public PotentialOrder1 {

            public:
                CubePotential(const Kernel * const kernel);

                virtual CubePotential *replicate() const override = 0;

                virtual void configureForType(const unsigned int &type) override;

                const Vec3 &getOrigin() const;
                void setOrigin(const Vec3 &origin);
                const Vec3 &getExtent() const;
                void setExtent(const Vec3 &extent);
                double getForceConstant() const;
                void setForceConstant(double forceConstant);
                bool isConsiderParticleRadius() const;
                void setConsiderParticleRadius(bool considerParticleRadius);
                double getParticleRadius() const;

            protected:
                const Kernel * const kernel;
                Vec3 origin;
                Vec3 extent;
                Vec3 min {}, max {};
                double particleRadius;
                double forceConstant;
                bool considerParticleRadius;
            };

            namespace _internal {
                template <> struct PotentialName<CubePotential> {static const std::string value;};
            }
        }
    }
}

#endif //READDY_MAIN_POTENTIALSORDER1_H
