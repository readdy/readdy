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
 * @file Reactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GPL-3
 */

#include <readdy/model/observables/Reactions.h>
#include <readdy/model/IOUtils.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>

namespace readdy {
namespace model {
namespace observables {

struct Reactions::Impl {
    using reactions_record_t = readdy::model::reactions::ReactionRecord;
    using reactions_writer_t = h5rd::VLENDataSet;
    std::unique_ptr<reactions_writer_t> writer;
    std::unique_ptr<util::TimeSeriesWriter> time;
    std::unique_ptr<h5rd::Group> group;
    std::unique_ptr<util::CompoundH5Types> h5types;
};

Reactions::Reactions(Kernel *const kernel, stride_type stride)
        : super(kernel, stride), pimpl(std::make_unique<Impl>()) {
}

void Reactions::flush() {
    if (pimpl->writer) pimpl->writer->flush();
    if (pimpl->time) pimpl->time->flush();
}

void Reactions::initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) {
    result.clear();
    h5rd::dimensions fs = {flushStride};
    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
    pimpl->h5types = std::make_unique<util::CompoundH5Types>(util::getReactionRecordTypes(file.ref()));
    pimpl->group = std::make_unique<h5rd::Group>(
            file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName));
    pimpl->writer = pimpl->group->createVLENDataSet("records", fs, dims, std::get<0>(*pimpl->h5types),
                                                    std::get<1>(*pimpl->h5types));
    pimpl->time = std::make_unique<util::TimeSeriesWriter>(*pimpl->group, flushStride);
}

void Reactions::append() {
    pimpl->writer->append({1}, &result);
    pimpl->time->append(t_current);
}

void Reactions::initialize(Kernel *const kernel) {
    if (!kernel->context().recordReactionsWithPositions()) {
        log::warn("The \"Reactions\"-observable set context.recordReactionsWithPositions() to true. "
                          "If this is undesired, the observable should not be registered.");
        kernel->context().recordReactionsWithPositions() = true;
    }
}

std::string Reactions::type() const {
    return "Reactions";
}

Reactions::~Reactions() = default;

}
}
}
