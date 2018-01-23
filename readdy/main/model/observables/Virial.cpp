/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * @file Virial.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 1/17/18
 */

#include <readdy/model/observables/Virial.h>

#include <readdy/model/Kernel.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>
#include <readdy/model/observables/io/Types.h>

namespace readdy {
namespace model {
namespace observables {

struct Virial::Impl {
    std::unique_ptr<h5rd::DataSet> ds{nullptr};
    std::unique_ptr<util::TimeSeriesWriter> time{nullptr};
    io::BloscFilter bloscFilter{};
};

Virial::Virial(Kernel *kernel, stride_type stride) : super(kernel, stride), pimpl(std::make_unique<Impl>()) {}

void Virial::flush() {
    if (pimpl->ds) pimpl->ds->flush();
    if (pimpl->time) pimpl->time->flush();
}

void Virial::initialize(Kernel *const kernel) {
    kernel->context().recordVirial() = true;
}

void Virial::initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) {
    h5rd::dimensions fs = {flushStride, Matrix33::n(), Matrix33::m()};
    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS, Matrix33::n(), Matrix33::m()};
    auto group = file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName);
    pimpl->ds = group.createDataSet<readdy::scalar>("data", fs, dims, {&pimpl->bloscFilter});
    pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
}

void Virial::append() {
    pimpl->ds->append({1, Matrix33::n(), Matrix33::m()}, result.data().data());
    pimpl->time->append(t_current);
}

std::string Virial::type() const {
    return "Virial";
}

Virial::~Virial() = default;


}
}
}
