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
 * @file PositionsObservable.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/observables/Positions.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>

namespace readdy {
namespace model {
namespace observables {

struct Positions::Impl {
    using writer_t = h5rd::VLENDataSet;
    std::unique_ptr<writer_t> writer;
    std::unique_ptr<util::TimeSeriesWriter> time;
    std::unique_ptr<util::CompoundH5Types> h5types;
};

Positions::Positions(Kernel *const kernel, stride_type stride,
                     std::vector<std::string> typesToCount) :
        Positions(kernel, stride,
                  _internal::util::transformTypes2(typesToCount, kernel->context())) {}

Positions::Positions(Kernel *const kernel, stride_type stride,
                     std::vector<ParticleTypeId> typesToCount) :
        Observable(kernel, stride), typesToCount(std::move(typesToCount)), pimpl(std::make_unique<Impl>()) {}

void Positions::append() {
    std::vector<Vec3> podVec(result.begin(), result.end());
    pimpl->writer->append({1}, &podVec);
    pimpl->time->append(t_current);
}

Positions::Positions(Kernel *const kernel, stride_type stride) : Observable(kernel, stride) {}

void Positions::initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) {
    pimpl->h5types = std::make_unique<util::CompoundH5Types>(util::getVec3Types(file.ref()));
    h5rd::dimensions fs = {flushStride};
    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
    auto group = file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName);
    pimpl->writer = group.createVLENDataSet("data", fs, dims, std::get<0>(*pimpl->h5types),
                                            std::get<1>(*pimpl->h5types));
    pimpl->time = std::make_unique<util::TimeSeriesWriter>(group, flushStride);
}

void Positions::flush() {
    if (pimpl->writer) pimpl->writer->flush();
    if (pimpl->time) pimpl->time->flush();
}

std::string Positions::type() const {
    return "Positions";
}

Positions::~Positions() = default;
}
}
}