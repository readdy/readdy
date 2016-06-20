/**
 * << detailed description >>
 *
 * @file P1Cube.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder1.h>

namespace readdy {

    namespace kernel {
        namespace singlecpu {
            namespace potentials {
                CubePotential::CubePotential(const readdy::kernel::singlecpu::SingleCPUKernel *const Kernel) : readdy::model::potentials::CubePotential<SingleCPUKernel>(Kernel) {

                }

                double CubePotential::calculateEnergy(const readdy::model::Vec3 &position) {

                    auto r = particleRadius;
                    if (!isConsiderParticleRadius()) r = 0;

                    double energy = 0;

                    for (auto i = 0; i < 3; ++i) {
                        if (position[i] - r < min[i] || position[i] + r > max[i]) {
                            if (position[i] - r < min[i]) {
                                energy += forceConstant * (position[i] - r - min[i]) * (position[i] - r - min[i]);
                            } else {
                                energy += forceConstant * (position[i] + r - max[i]) * (position[i] + r - max[i]);
                            }
                        }
                    }

                    return energy;
                }

                void CubePotential::calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) {

                    auto r = particleRadius;
                    if (!isConsiderParticleRadius()) r = 0;

                    for (auto i = 0; i < 3; i++) {
                        if (position[i] - r < min[i] || position[i] + r > max[i]) {
                            if (position[i] - r < min[i]) {
                                force[i] += forceConstant * (position[i] - r - min[i]);
                            } else {
                                force[i] += forceConstant * (position[i] + r - max[i]);
                            }
                        }
                    }

                }

                void CubePotential::calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &position) {

                    energy += calculateEnergy(position);
                    calculateForce(force, position);

                }

                potentials::CubePotential *CubePotential::replicate() const {
                    return new CubePotential(*this);
                }


            }
        }
    }

}


