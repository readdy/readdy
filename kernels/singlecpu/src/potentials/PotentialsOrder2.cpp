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

                double HarmonicRepulsion::calculateEnergy(const size_t &i, const size_t &j) {
                    const auto &&data = kernel->getKernelStateModelSingleCPU().getParticleData();
                    readdy::model::Vec3 positionI = *(data->begin_positions() + i);
                    readdy::model::Vec3 positionJ = *(data->begin_positions() + j);
                    double distance = (positionJ - positionI)*(positionJ - positionI);
                    if(distance < getSumOfParticleRadiiSquared()) {
                        distance = std::sqrt(distance);
                        distance -= getSumOfParticleRadii();
                        distance *= distance;
                        return distance * getForceConstant();
                    } else {
                        return 0;
                    }
                }

                void HarmonicRepulsion::calculateForce(readdy::model::Vec3 &force, const size_t &i, const size_t &j) {
                    const auto &&data = kernel->getKernelStateModelSingleCPU().getParticleData();
                    readdy::model::Vec3 positionI = *(data->begin_positions() + i);
                    readdy::model::Vec3 positionJ = *(data->begin_positions() + j);
                    readdy::model::Vec3 r_ij = positionJ - positionI;
                    double distance = r_ij * r_ij;
                    if (distance < getSumOfParticleRadiiSquared()) {
                        distance = std::sqrt(distance);
                        force = (2*getForceConstant() * (distance - getSumOfParticleRadii()))/distance * r_ij;
                    } else {
                        force *= 0;
                    }
                }

                void HarmonicRepulsion::calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const size_t &i, const size_t &j) {
                    const auto &&data = kernel->getKernelStateModelSingleCPU().getParticleData();
                    readdy::model::Vec3 positionI = *(data->begin_positions() + i);
                    readdy::model::Vec3 positionJ = *(data->begin_positions() + j);
                    readdy::model::Vec3 r_ij = positionJ - positionI;
                    double distance = r_ij * r_ij;
                    if (distance < getSumOfParticleRadiiSquared()) {
                        distance = std::sqrt(distance);
                        energy = getForceConstant() * std::pow(distance - getSumOfParticleRadii(), 2);
                        force = (2*getForceConstant() * (distance - getSumOfParticleRadii()))/distance * r_ij;
                    } else {
                        energy = 0;
                        force *= 0;
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