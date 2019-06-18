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
 * @file ReactionCounts.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright BSD-3
 */

#include <readdy/io/BloscFilter.h>
#include <readdy/model/observables/ReactionCounts.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>

namespace readdy::model::observables {

using data_set = h5rd::DataSet;
using ReactionsDataSets = std::unordered_map<ReactionId, std::unique_ptr<data_set>>;
using SpatialReactionsDataSets = std::unordered_map<ReactionId, std::unique_ptr<data_set>>;
using StructuralReactionsDataSets = std::unordered_map<ReactionId, std::unique_ptr<data_set>>;

struct ReactionCounts::Impl {
    std::unique_ptr<h5rd::Group> group;
    ReactionsDataSets dataSets;
    SpatialReactionsDataSets spatialReactionsDataSets;
    StructuralReactionsDataSets structuralReactionsDataSets;
    std::unique_ptr<util::TimeSeriesWriter> time;
    unsigned int flushStride = 0;
    bool firstWrite = true;
    std::function<void(std::unique_ptr<data_set> &)> flushFun = [](std::unique_ptr<data_set> &value) {
        if(value) value->flush();
    };
    io::BloscFilter bloscFilter {};
};

ReactionCounts::ReactionCounts(Kernel *const kernel, Stride stride)
        : Observable(kernel, stride), pimpl(std::make_unique<Impl>()) {
}

void ReactionCounts::flush() {
    readdy::util::collections::for_each_value_in_map(pimpl->dataSets, pimpl->flushFun);
    if (pimpl->time) pimpl->time->flush();
}

void ReactionCounts::initialize(Kernel *const kernel) {
    kernel->context().recordReactionCounts() = true;
}

void ReactionCounts::initializeDataSet(File &file, const std::string &dataSetName, Stride flushStride) {
    pimpl->firstWrite = true;
    pimpl->group = std::make_unique<h5rd::Group>(
            file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName));
    pimpl->flushStride = flushStride;
    pimpl->time = std::make_unique<util::TimeSeriesWriter>(*pimpl->group, flushStride);
}

void ReactionCounts::append() {
    const auto &[reactionCounts, spatialReactionCounts, structuralReactionCounts] = result;

    if (pimpl->firstWrite) {
        pimpl->firstWrite = false;

        {
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
        {
            const auto &topologyRegistry = kernel->context().topologyRegistry();
            auto spatialSubgroup = pimpl->group->createGroup("spatialCounts");
            for (const auto &entry : topologyRegistry.spatialReactionRegistry()) {
                for (const auto &spatialReaction : entry.second) {
                    h5rd::dimensions chunkSize = {pimpl->flushStride};
                    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
                    auto dset = spatialSubgroup.createDataSet<std::size_t>(std::to_string(spatialReaction.id()),
                            chunkSize, dims, {&pimpl->bloscFilter});
                    pimpl->spatialReactionsDataSets[spatialReaction.id()] = std::move(dset);
                }
            }
        }
        {
            const auto &topologyRegistry = kernel->context().topologyRegistry();
            auto structuralSubgroup = pimpl->group->createGroup("structuralCounts");
            for (const auto &entry : topologyRegistry.types()) {
                for (const auto &structuralReaction : entry.structuralReactions) {
                    h5rd::dimensions chunkSize = {pimpl->flushStride};
                    h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
                    auto dset = structuralSubgroup.createDataSet<std::size_t>(std::to_string(structuralReaction.id()),
                                                                           chunkSize, dims, {&pimpl->bloscFilter});
                    pimpl->structuralReactionsDataSets[structuralReaction.id()] = std::move(dset);
                }
            }
        }
    }
    for (const auto [id, count] : reactionCounts) {
        pimpl->dataSets.at(id)->append({1}, &count);
    }
    for(const auto [id, count] : spatialReactionCounts) {
        pimpl->spatialReactionsDataSets.at(id)->append({1}, &count);
    }
    for(const auto [id, count] : structuralReactionCounts) {
        pimpl->structuralReactionsDataSets.at(id)->append({1}, &count);
    }
    pimpl->time->append(t_current);
}

constexpr static auto& t = "ReactionCounts";

std::string_view ReactionCounts::type() const {
    return t;
}

ReactionCounts::~ReactionCounts() = default;

}
