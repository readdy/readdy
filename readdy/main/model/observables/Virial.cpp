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
