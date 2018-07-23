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
 * @file HistogramAlongAxis.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GPL-3
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