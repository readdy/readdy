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
 * @file PositionsObservable.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GPL-3
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