/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          *
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
 * @file Energy.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12/21/17
 */

#include "readdy/model/Kernel.h"
#include <readdy/model/observables/io/TimeSeriesWriter.h>
#include <readdy/model/observables/io/Types.h>
#include "readdy/model/observables/Energy.h"

namespace readdy {
namespace model {
namespace observables {

struct Energy::Impl {
    std::unique_ptr<h5rd::DataSet> ds{nullptr};
    std::unique_ptr<util::TimeSeriesWriter> time{nullptr};
    io::BloscFilter bloscFilter{};
};

Energy::Energy(Kernel *kernel, stride_type stride) : Observable(kernel, stride), pimpl(std::make_unique<Impl>()) {}

void Energy::flush() {
    if (pimpl->ds) pimpl->ds->flush();
    if (pimpl->time) pimpl->time->flush();
}

void Energy::initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) {
    h5rd::dimensions fs = {flushStride};
    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
    auto group = file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName);
    pimpl->ds = group.createDataSet<scalar>("data", fs, dims, {&pimpl->bloscFilter});
    pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
}

void Energy::append() {
    pimpl->ds->append({1}, &result);
    pimpl->time->append(t_current);
}

void Energy::evaluate() {
    result = kernel->stateModel().energy();
}

std::string Energy::type() const {
    return "Energy";
}

Energy::~Energy() = default;

}
}
}
