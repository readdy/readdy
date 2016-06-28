/**
 * << detailed description >>
 *
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#include <readdy/model/Observables.h>
#include <readdy/model/Kernel.h>

namespace readdy {
    namespace model {

        ObservableBase::~ObservableBase() {
            kernel->deconnectObservable(this);
        }

        void ParticlePositionObservable::evaluate() {
            ParticlePositionObservable::result = kernel->getKernelStateModel().getParticlePositions();
        }

        void TestCombinerObservable::evaluate() {
            std::vector<double> result;
            const auto &r1 = obs1->getResult();
            const auto &r2 = obs2->getResult();

            auto b1 = r1.begin();
            auto b2 = r2.begin();

            for (; b1 != r1.end();) {
                result.push_back((*b1) * (*b2));
                ++b1;
                ++b2;
            }

            TestCombinerObservable::result = result;
        }

        RadialDistributionObservable::RadialDistributionObservable(Kernel *const kernel, unsigned int stride, const std::vector<double> &binBorders, unsigned int typeCountFrom, unsigned int typeCountTo, double particleDensity)
                : Observable(kernel, stride), typeCountFrom(typeCountFrom), typeCountTo(typeCountTo), particleDensity(particleDensity) {
            setBinBorders(binBorders);
        }

        void RadialDistributionObservable::evaluate() {
            if(binBorders.size() > 1) {

                const auto particles = kernel->getKernelStateModel().getParticles();
                const auto n_particles = particles.size();
                {
                    const auto &distSquared = kernel->getKernelContext().getDistSquaredFun();
                    for (auto &&pFrom : particles) {
                        if (pFrom.getType() == typeCountFrom) {
                            for (auto &&pTo : particles) {
                                if (pTo.getType() == typeCountTo && pFrom.getId() != pTo.getId()) {
                                    const auto dist = sqrt(distSquared(pFrom.getPos(), pTo.getPos()));
                                    auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), dist);
                                    if (upperBound != binBorders.end()) {
                                        const auto binBordersIdx = upperBound - binBorders.begin();
                                        counts[binBordersIdx - 1]++;
                                    }
                                }
                            }
                        }
                    }
                }

                auto &radialDistribution = std::get<1>(result);
                {
                    const auto &binCenters = std::get<0>(result);
                    auto&& it_centers = binCenters.begin();
                    auto&& it_distribution = radialDistribution.begin();
                    for(auto&& it_counts = counts.begin(); it_counts != counts.end(); ++it_counts) {
                        const auto idx = it_centers - binCenters.begin();
                        const auto r = *it_centers;
                        const auto dr = binBorders[idx+1] - binBorders[idx];
                        *it_distribution = (*it_counts)/(4*M_PI*r*r*dr*n_particles*particleDensity);

                        ++it_distribution;
                        ++it_centers;
                    }
                }
            }
        }

        const std::vector<double> &RadialDistributionObservable::getBinBorders() const {
            return binBorders;
        }

        void RadialDistributionObservable::setBinBorders(const std::vector<double> &binBorders) {
            if(binBorders.size() > 1) {
                RadialDistributionObservable::binBorders = binBorders;
                auto nCenters = binBorders.size() - 1;
                result = std::make_pair(std::vector<double>(nCenters), std::vector<double>(nCenters));
                counts = std::vector<double>(nCenters);
                auto &binCenters = std::get<0>(result);
                auto it_begin = binBorders.begin();
                auto it_begin_next = it_begin + 1;
                size_t idx = 0;
                while (it_begin_next != binBorders.end()) {
                    binCenters[idx++] = (*it_begin + *it_begin_next) / 2;
                    ++it_begin;
                    ++it_begin_next;
                }
            } else {
                BOOST_LOG_TRIVIAL(warning) << "Argument bin borders' size should be at least two to make sense.";
            }

        }

    }


}




