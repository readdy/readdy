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
 * @file Forces.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/observables/Forces.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>

namespace readdy {
namespace model {
namespace observables {

struct Forces::Impl {
    using data_set_t = h5rd::VLENDataSet;
    std::unique_ptr<data_set_t> dataSet;
    std::unique_ptr<util::TimeSeriesWriter> timeSeries;
    std::unique_ptr<util::CompoundH5Types> h5types;
};

Forces::Forces(Kernel *const kernel, stride_type stride, std::vector<std::string> typesToCount)
        : Forces(kernel, stride,
                 _internal::util::transformTypes2(typesToCount, kernel->context())) {}

Forces::Forces(Kernel *const kernel, stride_type stride, const std::vector<ParticleTypeId> &typesToCount)
        : Observable(kernel, stride), typesToCount(typesToCount), pimpl(std::make_unique<Impl>()) {}

Forces::Forces(Kernel *const kernel, stride_type stride) : Observable(kernel, stride), typesToCount({}),
                                                            pimpl(std::make_unique<Impl>()) {}

void Forces::flush() {
    if (pimpl->dataSet) pimpl->dataSet->flush();
    if (pimpl->timeSeries) pimpl->timeSeries->flush();
}

void Forces::initializeDataSet(File &file, const std::string &dataSetName, unsigned int flushStride) {
    pimpl->h5types = std::make_unique<util::CompoundH5Types>(util::getVec3Types(file.parentFile()));
    h5rd::dimensions fs = {flushStride};
    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
    auto group = file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName);
    auto dataSet = group.createVLENDataSet("data", fs, dims,
                                           std::get<0>(*pimpl->h5types), std::get<1>(*pimpl->h5types));
    pimpl->dataSet = std::move(dataSet);
    pimpl->timeSeries = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
}

void Forces::append() {
    pimpl->dataSet->append({1}, &result);
    pimpl->timeSeries->append(t_current);
}

std::string Forces::type() const {
    return "Forces";
}

Forces::~Forces() = default;

}
}
}
