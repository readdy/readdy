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
 * @file RadialDistributionObservable.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/observables/RadialDistribution.h>

#include <utility>

#include <readdy/model/Kernel.h>
#include <readdy/common/numeric.h>
#include <readdy/model/observables/io/Types.h>

#include <readdy/model/observables/io/TimeSeriesWriter.h>
#include <readdy/common/boundary_condition_operations.h>

namespace readdy {
namespace model {
namespace observables {

struct RadialDistribution::Impl {
    using writer_t = h5rd::DataSet;
    std::unique_ptr<writer_t> writerRadialDistribution;
    std::unique_ptr<util::TimeSeriesWriter> time;
};

RadialDistribution::RadialDistribution(Kernel *const kernel, stride_type stride,
                                       std::vector<scalar> binBorders, std::vector<ParticleTypeId> typeCountFrom,
                                       std::vector<ParticleTypeId> typeCountTo, scalar particleToDensity)
        : Observable(kernel, stride), typeCountFrom(std::move(typeCountFrom)), typeCountTo(std::move(typeCountTo)),
          particleToDensity(particleToDensity), pimpl(std::make_unique<Impl>()) {
    setBinBorders(binBorders);
}

void RadialDistribution::evaluate() {
    if (binBorders.size() > 1) {
        std::fill(counts.begin(), counts.end(), 0);
        const auto particles = kernel->stateModel().getParticles();
        auto isInCollection = [](const readdy::model::Particle &p, const auto &collection) {
            return std::find(collection.begin(), collection.end(), p.getType()) != collection.end();
        };
        const auto nFromParticles = std::count_if(particles.begin(), particles.end(),
                                                  [this, isInCollection](const readdy::model::Particle &p) {
                                                      return isInCollection(p, typeCountFrom);
                                                  });
        {
            const auto pbc = kernel->context().periodicBoundaryConditions().data();
            const auto box = kernel->context().boxSize().data();
            for (auto &&pFrom : particles) {
                if (isInCollection(pFrom, typeCountFrom)) {
                    for (auto &&pTo : particles) {
                        if (isInCollection(pTo, typeCountTo) && pFrom.getId() != pTo.getId()) {
                            const auto dist = sqrt(bcs::distSquared(pFrom.getPos(), pTo.getPos(), box, pbc));
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
                        (c_::four / c_::three * readdy::util::numeric::pi<scalar>() * (std::pow(upperRadius, c_::three)
                                                                       - std::pow(lowerRadius, c_::three))
                         * nFromParticles * particleToDensity);
                ++it_distribution;
                ++it_centers;
            }
        }
    }
}

const std::vector<scalar> &RadialDistribution::getBinBorders() const {
    return binBorders;
}

void RadialDistribution::setBinBorders(const std::vector<scalar> &binBorders) {
    if (binBorders.size() > 1) {
        RadialDistribution::binBorders = binBorders;
        auto nCenters = binBorders.size() - 1;
        result = std::make_pair(std::vector<scalar>(nCenters), std::vector<scalar>(nCenters));
        counts = std::vector<scalar>(nCenters);
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
        log::warn("Argument bin borders' size should be at least two to make sense.");
    }
}

RadialDistribution::RadialDistribution(Kernel *const kernel, unsigned int stride,
                                       const std::vector<scalar> &binBorders,
                                       const std::vector<std::string> &typeCountFrom,
                                       const std::vector<std::string> &typeCountTo, scalar particleToDensity)
        : RadialDistribution(kernel, stride, binBorders,
                             _internal::util::transformTypes2(typeCountFrom, kernel->context()),
                             _internal::util::transformTypes2(typeCountTo, kernel->context()),
                             particleToDensity
) {}

void RadialDistribution::initializeDataSet(File &file, const std::string &dataSetName, unsigned int flushStride) {
    auto &centers = std::get<0>(result);
    h5rd::dimensions fs = {flushStride, centers.size()};
    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS, centers.size()};
    const auto path = std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName;
    auto group = file.createGroup(path);
    log::debug("created group with path {}", path);
    group.write("bin_centers", centers);
    pimpl->writerRadialDistribution = group.createDataSet<scalar>("distribution", fs, dims, {&bloscFilter});
    pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
}

void RadialDistribution::append() {
    auto &dist = std::get<1>(result);
    pimpl->writerRadialDistribution->append({1, dist.size()}, dist.data());
    pimpl->time->append(t_current);
}

void RadialDistribution::flush() {
    if (pimpl->writerRadialDistribution) pimpl->writerRadialDistribution->flush();
    if (pimpl->time) pimpl->time->flush();
}

std::string RadialDistribution::type() const {
    return "RadialDistribution";
}

RadialDistribution::~RadialDistribution() = default;
}
}
}