/**
 * << detailed description >>
 *
 * @file PotentialsOrder2.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/kernel/singlecpu/potentials/PotentialsOrder2.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace potentials {

                double SingleCPUHarmonicRepulsion::calculateEnergy(const vec_t &x_ij) {
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

                void SingleCPUHarmonicRepulsion::calculateForce(vec_t &force, const vec_t &x_ij) {
                    auto squared = x_ij * x_ij;
                    if (squared < getSumOfParticleRadiiSquared() && squared > 0) {
                        squared = std::sqrt(squared);
                        force = (2 * getForceConstant() * (squared - getSumOfParticleRadii())) / squared * x_ij;
                    } else {
                        force = {0,0,0};
                    }
                }

                void SingleCPUHarmonicRepulsion::calculateForceAndEnergy(vec_t &force, double &energy, const vec_t &x_ij) {
                    auto squared = x_ij * x_ij;
                    if (squared < getSumOfParticleRadiiSquared() && squared > 0) {
                        squared = std::sqrt(squared);
                        energy += getForceConstant() * std::pow(squared - getSumOfParticleRadii(), 2);
                        force = (2 * getForceConstant() * (squared - getSumOfParticleRadii())) / squared * x_ij;
                    } else {
                        force = {0,0,0};
                    }
                }

                SingleCPUHarmonicRepulsion::SingleCPUHarmonicRepulsion(const SingleCPUKernel *const kernel) : readdy::model::potentials::HarmonicRepulsion(kernel) { }

                potentials::SingleCPUHarmonicRepulsion *SingleCPUHarmonicRepulsion::replicate() const {
                    return new SingleCPUHarmonicRepulsion(*this);
                }

                double SingleCPUHarmonicRepulsion::getCutoffRadius() const {
                    return sumOfParticleRadii;
                }
                

                SingleCPUWeakInteractionPiecewiseHarmonic::SingleCPUWeakInteractionPiecewiseHarmonic(const readdy::model::Kernel *const kernel)
                        : readdy::model::potentials::WeakInteractionPiecewiseHarmonic(kernel) { }

                potentials::SingleCPUWeakInteractionPiecewiseHarmonic *SingleCPUWeakInteractionPiecewiseHarmonic::replicate() const {
                    return new SingleCPUWeakInteractionPiecewiseHarmonic(*this);
                }

                double SingleCPUWeakInteractionPiecewiseHarmonic::calculateEnergy(const readdy::model::Vec3 &x_ij) {
                    const auto dist = sqrt(x_ij * x_ij);
                    const auto len_part2 = noInteractionDistance - desiredParticleDistance;
                    if (dist < desiredParticleDistance) {
                        // repulsive as we are closer than the desired distance
                        return .5 * forceConstant * (dist - desiredParticleDistance) * (dist - desiredParticleDistance) - depthAtDesiredDistance;
                    } else {
                        // attractive as we are further (but not too far) apart than the desired distance
                        if (dist < desiredParticleDistance + .5 * len_part2) {
                            return .5 * depthAtDesiredDistance * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) * (dist - desiredParticleDistance) * (dist - desiredParticleDistance) -
                                   depthAtDesiredDistance;
                        } else {
                            // if we are not too far apart but still further than in the previous case, attractive
                            if (dist < noInteractionDistance) {
                                return -0.5 * depthAtDesiredDistance * (1 / (0.5 * len_part2)) * (1 / (0.5 * len_part2)) * (dist - noInteractionDistance) * (dist - noInteractionDistance);
                            }
                        }
                    }
                    return 0;
                }

                void SingleCPUWeakInteractionPiecewiseHarmonic::calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &x_ij) {
                    const auto dist = sqrt(x_ij * x_ij);
                    const auto len_part2 = noInteractionDistance - desiredParticleDistance;
                    double factor = 0;
                    if (dist < desiredParticleDistance) {
                        // repulsive as we are closer than the desired distance
                        factor = -1 * forceConstant * (desiredParticleDistance - dist);
                    } else {
                        // attractive as we are further (but not too far) apart than the desired distance
                        if (dist < desiredParticleDistance + .5 * len_part2) {
                            factor = -1 * depthAtDesiredDistance * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) * (desiredParticleDistance - dist);
                        } else {
                            // if we are not too far apart but still further than in the previous case, attractive
                            if (dist < noInteractionDistance) {
                                factor = depthAtDesiredDistance * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) * (noInteractionDistance - dist);
                            }
                        }
                    }
                    if(dist > 0 && factor != 0) {
                        force = factor * x_ij / dist;
                    } else {
                        force = {0,0,0};
                    }
                }

                void SingleCPUWeakInteractionPiecewiseHarmonic::calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &x_ij) {
                    energy += calculateEnergy(x_ij);
                    calculateForce(force, x_ij);
                }

                double SingleCPUWeakInteractionPiecewiseHarmonic::getCutoffRadius() const {
                    return noInteractionDistance;
                }


            }
        }
    }
}

