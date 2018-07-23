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
 * @file Forces.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GPL-3
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
