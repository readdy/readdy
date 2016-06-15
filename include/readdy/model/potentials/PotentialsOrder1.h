/**
 * << detailed description >>
 *
 * @file PotentialsOrder1.h
 * @brief << brief description >>
 * @author clonker
 * @date 15.06.16
 */

#include <readdy/model/KernelContext.h>
#include "PotentialOrder1.h"

#ifndef READDY_MAIN_POTENTIALSORDER1_H
#define READDY_MAIN_POTENTIALSORDER1_H

namespace readdy {
    namespace model {
        namespace potentials {
            template<typename KernelType>
            class CubePotential : public PotentialOrder1 {

            public:
                CubePotential(const KernelType * const kernel) :
                        PotentialOrder1(getPotentialName<CubePotential<KernelType>>()),
                        kernel(kernel)
                { }
                ~CubePotential() {
                    readdy::model::KernelContext ctx = kernel->getKernelContext();
                    ctx.deregisterPotential(getId());
                }

                virtual void configureForType(const unsigned int &type) override {
                    readdy::model::KernelContext ctx = kernel->getKernelContext();
                    particleRadius = ctx.getParticleRadius(type);
                    for (auto i = 0; i < 3; i++) {
                        if (origin[i] < origin[i] + extent[i]) {
                            min[i] = origin[i];
                            max[i] = origin[i] + extent[i];
                        } else {
                            min[i] = origin[i] + extent[i];
                            max[i] = origin[i];
                        }
                    }
                }


                const Vec3 &getOrigin() const {
                    return origin;
                }

                void setOrigin(const Vec3 &origin) {
                    CubePotential::origin = origin;
                }

                const Vec3 &getExtent() const {
                    return extent;
                }

                void setExtent(const Vec3 &extent) {
                    CubePotential::extent = extent;
                }

                double getForceConstant() const {
                    return forceConstant;
                }

                void setForceConstant(double forceConstant) {
                    CubePotential::forceConstant = forceConstant;
                }


                bool isConsiderParticleRadius() const {
                    return considerParticleRadius;
                }

                void setConsiderParticleRadius(bool considerParticleRadius) {
                    CubePotential::considerParticleRadius = considerParticleRadius;
                }

            protected:
                const KernelType * const kernel;
                Vec3 origin;
                Vec3 extent;
                Vec3 min {}, max {};
                double particleRadius;
                double forceConstant;
                bool considerParticleRadius;
            };

            namespace _internal {
                template<typename KernelType>
                struct PotentialName<CubePotential<KernelType>> {static const std::string value;};
            }
        }
    }
}

#include "_impl/PotentialsOrder1Impl.h"

#endif //READDY_MAIN_POTENTIALSORDER1_H
