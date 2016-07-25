/**
 * << detailed description >>
 *
 * @file CPUPotentials.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#include <readdy/kernel/cpu/potentials/CPUPotentials.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace potentials {
                CPUCubePotential::CPUCubePotential(const readdy::model::Kernel *const kernel)
                        : readdy::model::potentials::CubePotential(kernel) { }

                double CPUCubePotential::calculateEnergy(const readdy::model::Vec3 &position) {
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

                void CPUCubePotential::calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &position) {
                    auto r = particleRadius;
                    if (!isConsiderParticleRadius()) r = 0;
                    for (auto i = 0; i < 3; i++) {
                        if (position[i] - r < min[i] || position[i] + r > max[i]) {
                            if (position[i] - r < min[i]) {
                                force[i] += -1 * forceConstant * (position[i] - r - min[i]);
                            } else {
                                force[i] += -1 * forceConstant * (position[i] + r - max[i]);
                            }
                        }
                    }
                }

                void CPUCubePotential::calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy,
                                                               const readdy::model::Vec3 &position) {
                    energy += calculateEnergy(position);
                    calculateForce(force, position);
                }

                CPUCubePotential *CPUCubePotential::replicate() const {
                    return new CPUCubePotential(*this);
                }


            }
        }
    }
}

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace potentials {
                CPUHarmonicRepulsion::CPUHarmonicRepulsion(const readdy::model::Kernel *const kernel)
                        : HarmonicRepulsion(kernel) { }

                double CPUHarmonicRepulsion::calculateEnergy(const model::Vec3 &x_ij) {
                    auto distanceSquared = x_ij * x_ij;
                    if (distanceSquared < getSumOfParticleRadiiSquared()) {
                        distanceSquared = std::sqrt(distanceSquared);
                        distanceSquared -= getSumOfParticleRadii();
                        distanceSquared *= distanceSquared;
                        return distanceSquared * getForceConstant();
                    } else {
                        return 0;
                    }
                }

                void CPUHarmonicRepulsion::calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &x_ij) {
                    auto squared = x_ij * x_ij;
                    if (squared < getSumOfParticleRadiiSquared() && squared > 0) {
                        squared = std::sqrt(squared);
                        force = (2 * getForceConstant() * (squared - getSumOfParticleRadii())) / squared * x_ij;
                    } else {
                        force = {0, 0, 0};
                    }
                }

                void CPUHarmonicRepulsion::calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy,
                                                                   const readdy::model::Vec3 &x_ij) {
                    auto squared = x_ij * x_ij;
                    if (squared < getSumOfParticleRadiiSquared() && squared > 0) {
                        squared = std::sqrt(squared);
                        energy += getForceConstant() * std::pow(squared - getSumOfParticleRadii(), 2);
                        force = (2 * getForceConstant() * (squared - getSumOfParticleRadii())) / squared * x_ij;
                    } else {
                        force = {0, 0, 0};
                    }
                }

                CPUHarmonicRepulsion *CPUHarmonicRepulsion::replicate() const {
                    return new CPUHarmonicRepulsion(*this);
                }

                double CPUHarmonicRepulsion::getCutoffRadius() const {
                    return sumOfParticleRadii;
                }


            }
        }
    }
}

