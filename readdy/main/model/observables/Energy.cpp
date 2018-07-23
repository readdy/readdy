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
