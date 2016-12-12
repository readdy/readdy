/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file SingleCPUObservables.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */


#ifndef READDY_MAIN_SINGLECPUOBSERVABLES_H
#define READDY_MAIN_SINGLECPUOBSERVABLES_H

#include <readdy/model/observables/Observables.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>


namespace readdy {
namespace kernel {
namespace scpu {
namespace observables {

class SCPUPositions : public readdy::model::observables::Positions {
public:
    SCPUPositions(SCPUKernel *const kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {}) :
            readdy::model::observables::Positions(kernel, stride, typesToCount), kernel(kernel) {}

    virtual void evaluate() override {
        result.clear();
        const auto &pd = kernel->getKernelStateModel().getParticleData();
        auto positionsIt = pd->cbegin_positions();
        if (typesToCount.empty()) {
            // get all particles' positions
            result.reserve(pd->size());
            std::copy(positionsIt, pd->cend_positions(), std::back_inserter(result));
        } else {
            // only get positions of typesToCount
            auto typesIt = pd->cbegin_types();
            while (positionsIt != pd->cend_positions()) {
                if (std::find(typesToCount.begin(), typesToCount.end(), *typesIt) != typesToCount.end()) {
                    result.push_back(*positionsIt);
                }
                ++positionsIt;
                ++typesIt;
            }
        }
    }

protected:
    SCPUKernel *const kernel;
};

class SCPUParticles : public readdy::model::observables::Particles {
public:
    SCPUParticles(SCPUKernel *const kernel, unsigned int stride)
            : readdy::model::observables::Particles(kernel, stride), kernel(kernel) {};

    virtual void evaluate() override {
        auto &resultTypes = std::get<0>(result);
        auto &resultIds = std::get<1>(result);
        auto &resultPositions = std::get<2>(result);
        resultTypes.clear();
        resultIds.clear();
        resultPositions.clear();
        const auto &particleData = kernel->getKernelStateModel().getParticleData();
        auto typesIt = particleData->cbegin_types();
        auto idsIt = particleData->cbegin_ids();
        auto positionsIt = particleData->cbegin_positions();
        resultTypes.reserve(particleData->size());
        resultIds.reserve(particleData->size());
        resultPositions.reserve(particleData->size());
        std::copy(typesIt, particleData->cend_types(), std::back_inserter(resultTypes));
        std::copy(idsIt, particleData->cend_ids(), std::back_inserter(resultIds));
        std::copy(positionsIt, particleData->cend_positions(), std::back_inserter(resultPositions));
    };

protected:
    SCPUKernel *const kernel;
};

class SCPUHistogramAlongAxis : public readdy::model::observables::HistogramAlongAxis {

public:
    SCPUHistogramAlongAxis(SCPUKernel *const kernel, unsigned int stride,
                           const std::vector<double> &binBorders,
                           const std::vector<std::string> &typesToCount,
                           unsigned int axis)
            : readdy::model::observables::HistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis),
              kernel(kernel) {
        size = result.size();
    }

    virtual void evaluate() override {
        std::fill(result.begin(), result.end(), 0);

        const auto &model = kernel->getKernelStateModel();
        const auto data = model.getParticleData();

        auto it_pos = data->begin_positions();
        auto it_types = data->begin_types();

        while (it_pos != data->end_positions()) {
            if (typesToCount.find(*it_types) != typesToCount.end()) {
                const auto &vec = *it_pos;
                auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), vec[axis]);
                if (upperBound != binBorders.end()) {
                    unsigned long binBordersIdx = static_cast<unsigned long>(upperBound - binBorders.begin());
                    if (binBordersIdx > 1) {
                        ++result[binBordersIdx - 1];
                    }
                }
            }
            ++it_pos;
            ++it_types;
        }
    }

protected:
    SCPUKernel *const kernel;
    size_t size;
};

class SCPUNParticles : public readdy::model::observables::NParticles {
public:
    SCPUNParticles(SCPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {}) :
            readdy::model::observables::NParticles(kernel, stride, typesToCount),
            singleCPUKernel(kernel) {}

    virtual void evaluate() override {
        std::vector<unsigned long> resultVec = {};
        if (typesToCount.empty()) {
            resultVec.push_back(singleCPUKernel->getKernelStateModel().getParticleData()->size());
        } else {
            resultVec.resize(typesToCount.size());
            const auto &pd = singleCPUKernel->getKernelStateModel().getParticleData();
            auto typesIt = pd->cbegin_types();
            while (typesIt != pd->cend_types()) {
                unsigned int idx = 0;
                for (const auto t : typesToCount) {
                    if (*typesIt == t) {
                        resultVec[idx]++;
                        break;
                    }
                    ++idx;
                }
                ++typesIt;
            }
        }
        result = resultVec;
    }

protected:
    SCPUKernel *const singleCPUKernel;
};

class SCPUForces : public readdy::model::observables::Forces {
public:
    SCPUForces(SCPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {}) :
            readdy::model::observables::Forces(kernel, stride, typesToCount),
            kernel(kernel) {}

    virtual ~SCPUForces() {}

    virtual void evaluate() override {
        result.clear();
        const auto &pd = kernel->getKernelStateModel().getParticleData();
        auto forcesIt = pd->cbegin_forces();

        if (typesToCount.empty()) {
            // get all particles' forces
            result.reserve(pd->size());
            std::copy(forcesIt, pd->cend_forces(), std::back_inserter(result));
        } else {
            // only get forces of typesToCount
            auto typesIt = pd->cbegin_types();
            while (forcesIt != pd->cend_forces()) {
                for (auto countedParticleType : typesToCount) {
                    if (*typesIt == countedParticleType) {
                        result.push_back(*forcesIt);
                        break;
                    }
                }
                ++forcesIt;
                ++typesIt;
            }
        }
    }


protected:
    SCPUKernel *const kernel;
};

template<typename kernel_t=readdy::kernel::scpu::SCPUKernel>
class SCPURadialDistribution : public readdy::model::observables::RadialDistribution {
public:
    SCPURadialDistribution(kernel_t *const kernel, unsigned int stride, std::vector<double> binBorders, std::vector<std::string> typeCountFrom,
                                 std::vector<std::string> typeCountTo, double particleToDensity) :
            readdy::model::observables::RadialDistribution(kernel, stride, binBorders, typeCountFrom,
                                                           typeCountTo, particleToDensity), kernel(kernel) {}

    virtual void evaluate() override {
        if (binBorders.size() > 1) {
            std::fill(counts.begin(), counts.end(), 0);
            const auto particles = kernel->getKernelStateModel().getParticles();
            auto isInCollection = [](const readdy::model::Particle &p, const std::vector<unsigned int> &collection) {
                return std::find(collection.begin(), collection.end(), p.getType()) != collection.end();
            };
            const auto nFromParticles = std::count_if(particles.begin(), particles.end(),
                                                      [this, isInCollection](const readdy::model::Particle &p) {
                                                          return isInCollection(p, typeCountFrom);
                                                      });
            {
                const auto &distSquared = kernel->getKernelContext().getDistSquaredFun();
                for (auto &&pFrom : particles) {
                    if (isInCollection(pFrom, typeCountFrom)) {
                        for (auto &&pTo : particles) {
                            if (isInCollection(pTo, typeCountTo) && pFrom.getId() != pTo.getId()) {
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
                auto &&it_centers = binCenters.begin();
                auto &&it_distribution = radialDistribution.begin();
                for (auto &&it_counts = counts.begin(); it_counts != counts.end(); ++it_counts) {
                    const auto idx = it_centers - binCenters.begin();
                    const auto lowerRadius = binBorders[idx];
                    const auto upperRadius = binBorders[idx + 1];
                    *it_distribution =
                            (*it_counts) /
                            (4 / 3 * M_PI * (std::pow(upperRadius, 3) - std::pow(lowerRadius, 3)) * nFromParticles * particleToDensity);
                    ++it_distribution;
                    ++it_centers;
                }
            }
        }
    }

protected:
    kernel_t *const kernel;
};

}
}
}
}
#endif //READDY_MAIN_SINGLECPUOBSERVABLES_H
