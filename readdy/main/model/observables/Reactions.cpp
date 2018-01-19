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
 * @file Reactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
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
