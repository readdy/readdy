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
                HarmonicRepulsion::HarmonicRepulsion(const SingleCPUKernel *const kernel)
                        : readdy::model::potentials::HarmonicRepulsion<SingleCPUKernel>(readdy::model::potentials::_internal::PotentialName<HarmonicRepulsion>::value, kernel) {
                }

                double HarmonicRepulsion::calculateEnergy(const vec_t &x_i, const vec_t &x_j) {
                    auto distance = (x_j - x_i)*(x_j - x_i);
                    if(distance < getSumOfParticleRadiiSquared()) {
                        distance = std::sqrt(distance);
                        distance -= getSumOfParticleRadii();
                        distance *= distance;
                        return distance * getForceConstant();
                    } else {
                        return 0;
                    }
                }

                void HarmonicRepulsion::calculateForce(vec_t &force, const vec_t &x_i, const vec_t &x_j) {
                    const auto&& r_ij = x_j - x_i;
                    auto distance = r_ij * r_ij;
                    if (distance < getSumOfParticleRadiiSquared()) {
                        distance = std::sqrt(distance);
                        force += (2*getForceConstant() * (distance - getSumOfParticleRadii()))/distance * r_ij;
                    }
                }

                void HarmonicRepulsion::calculateForceAndEnergy(vec_t &force, double &energy, const vec_t &x_i, const vec_t &x_j) {
                    const auto&& r_ij = x_j - x_i;
                    auto distance = r_ij * r_ij;
                    if (distance < getSumOfParticleRadiiSquared()) {
                        distance = std::sqrt(distance);
                        energy += getForceConstant() * std::pow(distance - getSumOfParticleRadii(), 2);
                        force += (2*getForceConstant() * (distance - getSumOfParticleRadii()))/distance * r_ij;
                    }
                }
            }
        }
    }
    namespace model {
        namespace potentials {
            namespace _internal {
                const std::string PotentialName<readdy::kernel::singlecpu::potentials::HarmonicRepulsion>::value = "HarmonicRepulsion";
            }
        }
    }
}