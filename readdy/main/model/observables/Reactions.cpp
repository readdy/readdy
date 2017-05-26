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
    using reactions_writer_t = io::VLENDataSet;
    std::unique_ptr<reactions_writer_t> writer;
    std::unique_ptr<util::TimeSeriesWriter> time;
    std::unique_ptr<readdy::io::Group> group;
    bool firstWrite = true;
};

Reactions::Reactions(Kernel *const kernel, unsigned int stride)
        : super(kernel, stride), pimpl(std::make_unique<Impl>()) {
}

void Reactions::flush() {
    if (pimpl->writer) pimpl->writer->flush();
    if (pimpl->time) pimpl->time->flush();
}

void Reactions::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->writer) {
        std::vector<readdy::io::h5::dims_t> fs = {flushStride};
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        pimpl->group = std::make_unique<io::Group>(
                file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName));
        {
            auto dataSet = std::make_unique<Impl::reactions_writer_t>(
                    pimpl->group->createVLENDataSet("records", fs, dims, util::ReactionRecordPODMemoryType(),
                                                    util::ReactionRecordPODFileType()));
            pimpl->writer = std::move(dataSet);
        }
        pimpl->time = std::make_unique<util::TimeSeriesWriter>(*pimpl->group, flushStride);
    }
}

void Reactions::append() {
    pimpl->writer->append({1}, &result);
    if (pimpl->firstWrite) {
        pimpl->firstWrite = false;
        writeReactionInformation(*pimpl->group, kernel->getKernelContext());
    }
    pimpl->time->append(t_current);
}

void Reactions::initialize(Kernel *const kernel) {
    if (!kernel->getKernelContext().recordReactionsWithPositions()) {
        log::warn("The \"Reactions\"-observable set context.recordReactionsWithPositions() to true. "
                          "If this is undesired, the observable should not be registered.");
        kernel->getKernelContext().recordReactionsWithPositions() = true;
    }
}

Reactions::~Reactions() = default;

}
}
}