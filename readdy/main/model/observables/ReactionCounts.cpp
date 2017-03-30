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
    std::unique_ptr<data_set_t> ds_order1;
    std::unique_ptr<data_set_t> ds_order2;
    std::unique_ptr<util::TimeSeriesWriter> time;
    bool shouldWrite = false;
    unsigned int flushStride = 0;
    bool firstWrite = true;
};

ReactionCounts::ReactionCounts(Kernel *const kernel, unsigned int stride)
        : Observable(kernel, stride), pimpl(std::make_unique<Impl>()) {
    result = std::make_tuple(std::vector<std::size_t>(), std::vector<std::size_t>());
}

void ReactionCounts::flush() {
    if (pimpl->ds_order1) pimpl->ds_order1->flush();
    if (pimpl->ds_order2) pimpl->ds_order2->flush();
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
        const auto n_reactions_order1 = kernel->getKernelContext().reactions().n_order1();
        const auto n_reactions_order2 = kernel->getKernelContext().reactions().n_order2();
        std::get<0>(result).resize(n_reactions_order1);
        std::get<1>(result).resize(n_reactions_order2);
        if (pimpl->shouldWrite) {
            auto subgroup = pimpl->group->createGroup("counts");
            if(n_reactions_order1 > 0){
                std::vector<readdy::io::h5::dims_t> fs = {pimpl->flushStride, n_reactions_order1};
                std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS, n_reactions_order1};
                auto dataSetTypes = std::make_unique<Impl::data_set_t>("order1", subgroup, fs, dims);
                pimpl->ds_order1 = std::move(dataSetTypes);
            }
            if(n_reactions_order2 > 0){
                std::vector<readdy::io::h5::dims_t> fs = {pimpl->flushStride, n_reactions_order2};
                std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS, n_reactions_order2};
                auto dataSetTypes = std::make_unique<Impl::data_set_t>("order2", subgroup, fs, dims);
                pimpl->ds_order2 = std::move(dataSetTypes);
            }
            writeReactionInformation(*pimpl->group, kernel->getKernelContext());
        }
        pimpl->firstWrite = false;
    }
    if(pimpl->ds_order1){
        auto &order1Reactions = std::get<0>(result);
        pimpl->ds_order1->append({1, order1Reactions.size()}, order1Reactions.data());
    }
    if(pimpl->ds_order2){
        auto &order2Reactions = std::get<1>(result);
        pimpl->ds_order2->append({1, order2Reactions.size()}, order2Reactions.data());
    }
    pimpl->time->append(t_current);
}

ReactionCounts::~ReactionCounts() = default;

}
}
}