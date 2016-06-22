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
                    auto distanceSquared = x_ij*x_ij;
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
                    if (squared < getSumOfParticleRadiiSquared()) {
                        squared = std::sqrt(squared);
                        force += (2 * getForceConstant() * (squared - getSumOfParticleRadii())) / squared * x_ij;
                    }
                }

                void SingleCPUHarmonicRepulsion::calculateForceAndEnergy(vec_t &force, double &energy, const vec_t &x_ij) {
                    auto squared = x_ij * x_ij;
                    if (squared < getSumOfParticleRadiiSquared()) {
                        squared = std::sqrt(squared);
                        energy += getForceConstant() * std::pow(squared - getSumOfParticleRadii(), 2);
                        force += (2 * getForceConstant() * (squared - getSumOfParticleRadii())) / squared * x_ij;
                    }
                }

                SingleCPUHarmonicRepulsion::SingleCPUHarmonicRepulsion(const SingleCPUKernel *const kernel) : readdy::model::potentials::HarmonicRepulsion(kernel) { }

                potentials::SingleCPUHarmonicRepulsion *SingleCPUHarmonicRepulsion::replicate() const {
                    return new SingleCPUHarmonicRepulsion(*this);
                }


            }
        }
    }
}
