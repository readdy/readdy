/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file RadialDistributionObservable.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GPL-3
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
            return std::find(collection.begin(), collection.end(), p.type()) != collection.end();
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
                        if (isInCollection(pTo, typeCountTo) && pFrom.id() != pTo.id()) {
                            const auto dist = sqrt(bcs::distSquared(pFrom.pos(), pTo.pos(), box, pbc));
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