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
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 26.04.16
 */

#include <readdy/model/observables/Observables.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/_internal/Util.h>

namespace readdy {
namespace model {
namespace observables {

ParticlePosition::ParticlePosition(Kernel *const kernel, unsigned int stride,
                                   std::vector<std::string> typesToCount) :
        ParticlePosition(kernel, stride,
                         _internal::util::transformTypes2(typesToCount, kernel->getKernelContext())) {}

ParticlePosition::ParticlePosition(Kernel *const kernel, unsigned int stride,
                                   std::vector<unsigned int> typesToCount) :
        Observable(kernel, stride), typesToCount(typesToCount) {}

void TestCombiner::evaluate() {
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

    TestCombiner::result = result;
}

RadialDistribution::RadialDistribution(Kernel *const kernel, unsigned int stride,
                                       std::vector<double> binBorders, unsigned int typeCountFrom,
                                       unsigned int typeCountTo, double particleDensity)
        : Observable(kernel, stride), typeCountFrom(typeCountFrom), typeCountTo(typeCountTo),
          particleDensity(particleDensity) {
    setBinBorders(binBorders);
}

void RadialDistribution::evaluate() {
    if (binBorders.size() > 1) {
        std::fill(counts.begin(), counts.end(), 0);
        const auto particles = kernel->getKernelStateModel().getParticles();
        const auto n_from_particles = std::count_if(particles.begin(), particles.end(),
                                                    [this](const readdy::model::Particle &p) {
                                                        return p.getType() == typeCountFrom;
                                                    });
        {
            const auto &distSquared = kernel->getKernelContext().getDistSquaredFun();
            for (auto &&pFrom : particles) {
                if (pFrom.getType() == typeCountFrom) {
                    for (auto &&pTo : particles) {
                        if (pTo.getType() == typeCountTo && pFrom.getId() != pTo.getId()) {
                            const auto dist = std::sqrt(distSquared(pFrom.getPos(), pTo.getPos()));
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
                const auto r = *it_centers;
                const auto dr = binBorders[idx + 1] - binBorders[idx];
                *it_distribution = (*it_counts) / (4 * M_PI * r * r * dr * n_from_particles * particleDensity);

                ++it_distribution;
                ++it_centers;
            }
        }
    }
}

const std::vector<double> &RadialDistribution::getBinBorders() const {
    return binBorders;
}

void RadialDistribution::setBinBorders(const std::vector<double> &binBorders) {
    if (binBorders.size() > 1) {
        RadialDistribution::binBorders = binBorders;
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
        log::console()->warn("Argument bin borders' size should be at least two to make sense.");
    }

}

RadialDistribution::RadialDistribution(Kernel *const kernel, unsigned int stride,
                                       std::vector<double> binBorders,
                                       const std::string &typeCountFrom,
                                       const std::string &typeCountTo, double particleDensity)
        : RadialDistribution(
        kernel, stride, binBorders,
        kernel->getKernelContext().getParticleTypeID(typeCountFrom),
        kernel->getKernelContext().getParticleTypeID(typeCountTo),
        particleDensity
) {

}

CenterOfMass::CenterOfMass(readdy::model::Kernel *const kernel, unsigned int stride,
                           unsigned int particleType)
        : Observable(kernel, stride), particleTypes({particleType}) {
}

CenterOfMass::CenterOfMass(Kernel *const kernel, unsigned int stride,
                           const std::string &particleType)
        : CenterOfMass(kernel, stride, kernel->getKernelContext().getParticleTypeID(particleType)) {}

void CenterOfMass::evaluate() {
    readdy::model::Vec3 com{0, 0, 0};
    unsigned long n_particles = 0;
    for (auto &&p : kernel->getKernelStateModel().getParticles()) {
        if (particleTypes.find(p.getType()) != particleTypes.end()) {
            ++n_particles;
            com += p.getPos();
        }
    }
    com /= n_particles;
    result = com;
}

CenterOfMass::CenterOfMass(Kernel *const kernel, unsigned int stride,
                           const std::vector<unsigned int> &particleTypes)
        : Observable(kernel, stride), particleTypes(particleTypes.begin(), particleTypes.end()) {
}

CenterOfMass::CenterOfMass(Kernel *const kernel, unsigned int stride,
                           const std::vector<std::string> &particleType) : Observable(kernel,
                                                                                      stride),
                                                                           particleTypes() {
    for (auto &&pt : particleType) {
        particleTypes.emplace(kernel->getKernelContext().getParticleTypeID(pt));
    }

}

NParticles::NParticles(Kernel *const kernel, unsigned int stride,
                       std::vector<std::string> typesToCount)
        : NParticles(kernel, stride,
                     _internal::util::transformTypes2(typesToCount, kernel->getKernelContext())) {

}

NParticles::NParticles(Kernel *const kernel, unsigned int stride,
                       std::vector<unsigned int> typesToCount)
        : Observable(kernel, stride), typesToCount(typesToCount) {
}

Forces::Forces(Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount)
        : Forces(kernel, stride,
                 _internal::util::transformTypes2(typesToCount, kernel->getKernelContext())) {}

Forces::Forces(Kernel *const kernel, unsigned int stride, std::vector<unsigned int> typesToCount)
        : Observable(kernel, stride), typesToCount(typesToCount) {}

HistogramAlongAxis::HistogramAlongAxis(readdy::model::Kernel *const kernel, unsigned int stride,
                                       std::vector<double> binBorders,
                                       std::set<unsigned int> typesToCount, unsigned int axis)
        : Observable(kernel, stride), binBorders(binBorders), typesToCount(typesToCount), axis(axis) {
    auto nCenters = binBorders.size() - 1;
    result = std::vector<double>(nCenters);
}


HistogramAlongAxis::HistogramAlongAxis(Kernel *const kernel, unsigned int stride,
                                       std::vector<double> binBorders,
                                       std::vector<std::string> typesToCount,
                                       unsigned int axis)
        : HistogramAlongAxis(kernel, stride, binBorders,
                             _internal::util::transformTypes(typesToCount, kernel->getKernelContext()),
                             axis) {

}


}
}
}
