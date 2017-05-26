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
 * @file NParticles.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/observables/NParticles.h>
#include <readdy/model/Kernel.h>
#include <readdy/io/DataSet.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>

namespace readdy {
namespace model {
namespace observables {

NParticles::NParticles(Kernel *const kernel, unsigned int stride,
                       std::vector<std::string> typesToCount)
        : NParticles(kernel, stride,
                     _internal::util::transformTypes2(typesToCount, kernel->getKernelContext())) {

}

NParticles::NParticles(Kernel *const kernel, unsigned int stride,
                       std::vector<unsigned int> typesToCount)
        : Observable(kernel, stride), typesToCount(typesToCount), pimpl(std::make_unique<Impl>()) {
}

NParticles::~NParticles() = default;

struct NParticles::Impl {
    std::unique_ptr<io::DataSet> ds;
    std::unique_ptr<util::TimeSeriesWriter> time;
};

NParticles::NParticles(Kernel *const kernel, unsigned int stride)
        : Observable(kernel, stride), pimpl(std::make_unique<Impl>()) {}

void NParticles::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->ds) {
        const auto size = typesToCount.empty() ? 1 : typesToCount.size();
        std::vector<readdy::io::h5::dims_t> fs = {flushStride, size};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS, size};
        auto group = file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName);
        pimpl->ds = std::make_unique<io::DataSet>(group.createDataSet<std::size_t>("data", fs, dims));
        pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
    }
}

void NParticles::append() {
    pimpl->ds->append({1, result.size()}, result.data());
    pimpl->time->append(t_current);
}

void NParticles::flush() {
    if (pimpl->ds) pimpl->ds->flush();
    if(pimpl->time) pimpl->time->flush();
}
}
}
}