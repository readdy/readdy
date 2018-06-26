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
 * @file HistogramAlongAxis.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/observables/HistogramAlongAxis.h>

#include <readdy/model/Kernel.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>

namespace readdy {
namespace model {
namespace observables {


struct HistogramAlongAxis::Impl {
    std::unique_ptr<h5rd::DataSet> dataSet;
    std::unique_ptr<util::TimeSeriesWriter> time;
};

HistogramAlongAxis::HistogramAlongAxis(readdy::model::Kernel *const kernel, stride_type stride,
                                       std::vector<scalar> binBorders, std::set<ParticleTypeId> typesToCount,
                                       unsigned int axis)
        : Observable(kernel, stride), binBorders(binBorders), typesToCount(std::move(typesToCount)), axis(axis),
          pimpl(std::make_unique<Impl>()) {
    auto nCenters = binBorders.size() - 1;
    result = std::vector<scalar>(nCenters);
}


HistogramAlongAxis::HistogramAlongAxis(Kernel *const kernel, stride_type stride,
                                       std::vector<scalar> binBorders,
                                       std::vector<std::string> typesToCount,
                                       unsigned int axis)
        : HistogramAlongAxis(kernel, stride, std::move(binBorders),
                             _internal::util::transformTypes(typesToCount, kernel->context()),
                             axis) {

}

void HistogramAlongAxis::initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) {
    const auto size = result.size();
    h5rd::dimensions fs = {flushStride, size};
    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS, size};
    const auto path = std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName;
    auto group = file.createGroup(path);
    pimpl->dataSet = group.createDataSet<scalar>("data", fs, dims, {&bloscFilter});
    pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
}

void HistogramAlongAxis::append() {
    pimpl->dataSet->append({1, result.size()}, result.data());
    pimpl->time->append(t_current);
}

void HistogramAlongAxis::flush() {
    if (pimpl->dataSet) pimpl->dataSet->flush();
    if (pimpl->time) pimpl->time->flush();
}

std::string HistogramAlongAxis::type() const {
    return "HistogramAlongAxis";
}

HistogramAlongAxis::~HistogramAlongAxis() = default;

}
}
}