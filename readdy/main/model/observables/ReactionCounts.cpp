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
 * @file ReactionCounts.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/io/BloscFilter.h>
#include <readdy/model/observables/ReactionCounts.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>

namespace readdy {
namespace model {
namespace observables {

using data_set = h5rd::DataSet;
using data_set_map = std::unordered_map<reactions::Reaction::ReactionId, std::unique_ptr<data_set>>;

struct ReactionCounts::Impl {
    std::unique_ptr<h5rd::Group> group;
    data_set_map dataSets;
    std::unique_ptr<util::TimeSeriesWriter> time;
    unsigned int flushStride = 0;
    bool firstWrite = true;
    std::function<void(std::unique_ptr<data_set> &)> flushFun = [](std::unique_ptr<data_set> &value) {
        if(value) value->flush();
    };
    io::BloscFilter bloscFilter {};
};

ReactionCounts::ReactionCounts(Kernel *const kernel, stride_type stride)
        : Observable(kernel, stride), pimpl(std::make_unique<Impl>()) {
}

void ReactionCounts::flush() {
    readdy::util::collections::for_each_value_in_map(pimpl->dataSets, pimpl->flushFun);
    if (pimpl->time) pimpl->time->flush();
}

void ReactionCounts::initialize(Kernel *const kernel) {
    if (!kernel->context().recordReactionCounts()) {
        log::warn("The \"ReactionCounts\"-observable set context.recordReactionCounts() to true. "
                          "If this is undesired, the observable should not be registered.");
        kernel->context().recordReactionCounts() = true;
    }
}

void ReactionCounts::initializeDataSet(File &file, const std::string &dataSetName, stride_type flushStride) {
    pimpl->firstWrite = true;
    pimpl->group = std::make_unique<h5rd::Group>(
            file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName));
    pimpl->flushStride = flushStride;
    pimpl->time = std::make_unique<util::TimeSeriesWriter>(*pimpl->group, flushStride);
}

void ReactionCounts::append() {
    if (pimpl->firstWrite) {
        pimpl->firstWrite = false;
        const auto &reactionRegistry = kernel->context().reactions();
        auto subgroup = pimpl->group->createGroup("counts");
        for (const auto &reaction : reactionRegistry.order1Flat()) {
            h5rd::dimensions chunkSize = {pimpl->flushStride};
            h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
            auto dset = subgroup.createDataSet<std::size_t>(std::to_string(reaction->id()), chunkSize, dims,
                                                                 {&pimpl->bloscFilter});
            pimpl->dataSets[reaction->id()] = std::move(dset);
        }
        for (const auto &reaction : reactionRegistry.order2Flat()) {
            h5rd::dimensions chunkSize = {pimpl->flushStride};
            h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
            auto dset = subgroup.createDataSet<std::size_t>(std::to_string(reaction->id()), chunkSize, dims,
                                                            {&pimpl->bloscFilter});
            pimpl->dataSets[reaction->id()] = std::move(dset);
        }
    }
    for (const auto &reactionEntry : result) {
        auto copy = reactionEntry.second;
        pimpl->dataSets.at(reactionEntry.first)->append({1}, &copy);
    }
    pimpl->time->append(t_current);
}

std::string ReactionCounts::type() const {
    return "ReactionCounts";
}

ReactionCounts::~ReactionCounts() = default;

}
}
}
