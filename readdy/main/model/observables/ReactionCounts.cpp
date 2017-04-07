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

#include <readdy/model/observables/ReactionCounts.h>
#include <readdy/model/IOUtils.h>
#include <readdy/model/Kernel.h>
#include <readdy/io/DataSet.h>
#include <readdy/model/observables/io/Types.h>
#include <readdy/model/observables/io/TimeSeriesWriter.h>

namespace readdy {
namespace model {
namespace observables {


struct ReactionCounts::Impl {
    using data_set_t = io::DataSet<std::size_t, false>;
    std::unique_ptr<io::Group> group;
    std::unordered_map<readdy::particle_type_type, data_set_t> ds_order1;
    std::unordered_map<readdy::util::particle_type_pair, data_set_t,
            readdy::util::particle_type_pair_hasher, readdy::util::particle_type_pair_equal_to> ds_order2;
    std::unique_ptr<util::TimeSeriesWriter> time;
    bool shouldWrite = false;
    unsigned int flushStride = 0;
    bool firstWrite = true;
    std::function<void(data_set_t&)> flushFun = [](data_set_t &value){
        value.flush();
    };
};

ReactionCounts::ReactionCounts(Kernel *const kernel, unsigned int stride)
        : Observable(kernel, stride), pimpl(std::make_unique<Impl>()) {
}

void ReactionCounts::flush() {
    readdy::util::collections::for_each_value_in_map(pimpl->ds_order1, pimpl->flushFun);
    readdy::util::collections::for_each_value_in_map(pimpl->ds_order2, pimpl->flushFun);
    if (pimpl->time) pimpl->time->flush();
}

void ReactionCounts::initialize(Kernel *const kernel) {
    if (!kernel->getKernelContext().recordReactionCounts()) {
        log::warn("The \"ReactionCounts\"-observable set context.recordReactionCounts() to true. "
                          "If this is undesired, the observable should not be registered.");
        kernel->getKernelContext().recordReactionCounts() = true;
    }
}

void ReactionCounts::initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) {
    if (!pimpl->group) {
        pimpl->group = std::make_unique<io::Group>(
                file.createGroup(std::string(util::OBSERVABLES_GROUP_PATH) + "/" + dataSetName));
        pimpl->shouldWrite = true;
        pimpl->flushStride = flushStride;
        pimpl->time = std::make_unique<util::TimeSeriesWriter>(*pimpl->group, flushStride);
    }
}

void ReactionCounts::append() {
    if (pimpl->firstWrite) {
        const auto &ctx = kernel->getKernelContext();
        if (pimpl->shouldWrite) {
            auto subgroup = pimpl->group->createGroup("counts");
            auto order1Subgroup = subgroup.createGroup("order1");
            for (const auto &entry : ctx.particle_types().type_mapping()) {
                const auto &pType = entry.second;
                const auto numberOrder1Reactions = ctx.reactions().order1_by_type(pType).size();
                if (numberOrder1Reactions > 0) {
                    std::vector<readdy::io::h5::dims_t> chunkSize = {pimpl->flushStride, numberOrder1Reactions};
                    std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS, numberOrder1Reactions};
                    pimpl->ds_order1.emplace(std::piecewise_construct, std::forward_as_tuple(pType), std::forward_as_tuple(
                            ctx.particle_types().name_of(pType) + "[id=" + std::to_string(pType) + "]",
                            order1Subgroup, chunkSize, dims));
                }
            }
            auto order2Subgroup = subgroup.createGroup("order2");
            for (const auto &entry1 : ctx.particle_types().type_mapping()) {
                const auto &pType1 = entry1.second;
                for (const auto &entry2 : ctx.particle_types().type_mapping()) {
                    const auto &pType2 = entry2.second;
                    if (pType2 < pType1) continue;
                    const auto numberOrder2Reactions = ctx.reactions().order2_by_type(pType1, pType2).size();
                    if (numberOrder2Reactions > 0) {
                        std::vector<readdy::io::h5::dims_t> chunkSize = {pimpl->flushStride, numberOrder2Reactions};
                        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS, numberOrder2Reactions};
                        pimpl->ds_order2.emplace(std::piecewise_construct, std::forward_as_tuple(std::tie(pType1, pType2)), std::forward_as_tuple(
                                ctx.particle_types().name_of(pType1) + "[id=" + std::to_string(pType1) + "] + " +
                                ctx.particle_types().name_of(pType2) + "[id=" + std::to_string(pType2) + "]",
                                order2Subgroup, chunkSize, dims));
                    }
                }
            }
            writeReactionInformation(*pimpl->group, kernel->getKernelContext());
        }
        pimpl->firstWrite = false;
    }

    // actual writing of data
    const auto &countsOrder1 = std::get<0>(result);
    writeCountsToDataSets(countsOrder1, pimpl->ds_order1);

    const auto &countsOrder2 = std::get<1>(result);
    writeCountsToDataSets(countsOrder2, pimpl->ds_order2);

    pimpl->time->append(t_current);
}

ReactionCounts::~ReactionCounts() = default;

}
}
}