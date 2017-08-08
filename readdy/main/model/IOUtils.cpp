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
 * Write datasets that describe the context of the system that was used.
 * Since these datasets are meta-data, turn off compression, i.e. no blosc-filter is necessary for reading this.
 *
 * @file IOUtils.cpp
 * @brief Utils for writing to files that concern the model, e.g. writing reaction information.
 * @author clonker
 * @author chrisfroe
 * @date 10.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/IOUtils.h>
#include <iostream>

namespace readdy {
namespace model {
namespace ioutils {

void writeReactionInformation(readdy::io::Group &group, const KernelContext &context) {
    auto subgroup = group.createGroup("./registered_reactions");
    // order1
    const auto &order1_reactions = context.reactions().order1_flat();
    auto n_reactions = order1_reactions.size();
    if (n_reactions > 0) {
        std::vector<ReactionInfo> order1_info;
        for (const auto &r : order1_reactions) {
            const auto &reactions_current_type = context.reactions().order1_by_type(r->getEducts()[0]);
            auto it = std::find_if(reactions_current_type.begin(), reactions_current_type.end(),
                                   [&r](const reactions::Reaction<1> *x) { return x->getId() == r->getId(); });
            if (it != reactions_current_type.end()) {
                std::size_t index = static_cast<std::size_t>(it - reactions_current_type.begin());
                const std::array<particle_type_type, 2> educts = {r->getEducts()[0], 0};
                ReactionInfo info{r->getName().c_str(), index, r->getId(), r->getNEducts(), r->getNProducts(),
                                  r->getRate(), r->getEductDistance(),
                                  r->getProductDistance(), educts, r->getProducts()};
                order1_info.push_back(info);
            }
        }
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        std::vector<readdy::io::h5::dims_t> extent = {n_reactions};
        auto order1_reaction_dset = subgroup.createDataSet("order1_reactions", extent, dims, ReactionInfoMemoryType(),
                                                           ReactionInfoFileType(), io::DataSetCompression::none);
        order1_reaction_dset.append(extent, order1_info.data());
    }
    // order2
    const auto &order2_reactions = context.reactions().order2_flat();
    n_reactions = order2_reactions.size();
    if (n_reactions > 0) {
        std::vector<ReactionInfo> order2_info;
        for (const auto &r : order2_reactions) {
            const auto &reactions_current_type = context.reactions().order2_by_type(r->getEducts()[0],
                                                                                    r->getEducts()[1]);
            auto it = std::find_if(reactions_current_type.begin(), reactions_current_type.end(),
                                   [&r](const reactions::Reaction<2> *x) { return x->getId() == r->getId(); });
            if (it != reactions_current_type.end()) {
                std::size_t index = static_cast<std::size_t>(it - reactions_current_type.begin());
                ReactionInfo info{r->getName().c_str(), index, r->getId(), r->getNEducts(), r->getNProducts(),
                                  r->getRate(), r->getEductDistance(),
                                  r->getProductDistance(), r->getEducts(), r->getProducts()};
                order2_info.push_back(info);
            }
        }
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        std::vector<readdy::io::h5::dims_t> extent = {n_reactions};
        auto order2_reaction_dset = subgroup.createDataSet("order2_reactions", extent, dims, ReactionInfoMemoryType(),
                                                           ReactionInfoFileType(), io::DataSetCompression::none);
        order2_reaction_dset.append(extent, order2_info.data());
    }
}

void writeParticleTypeInformation(readdy::io::Group &group, const KernelContext &context) {
    const auto &types = context.particle_types().type_mapping();
    std::vector<ParticleTypeInfo> type_info_vec;
    for (const auto &p_type : types) {
        ParticleTypeInfo info{p_type.first.c_str(), p_type.second,
                              context.particle_types().diffusion_constant_of(p_type.first)};
        type_info_vec.push_back(info);
    }
    if (!type_info_vec.empty()) {
        std::vector<readdy::io::h5::dims_t> dims = {readdy::io::h5::UNLIMITED_DIMS};
        std::vector<readdy::io::h5::dims_t> extent = {type_info_vec.size()};
        auto dset = group.createDataSet("particle_types", extent, dims, ParticleTypeInfoMemoryType(),
                                        ParticleTypeInfoFileType(), io::DataSetCompression::none);
        dset.append(extent, type_info_vec.data());
    }
}

void writeSimulationSetup(io::Group &group, const KernelContext &context) {
    writeParticleTypeInformation(group, context);
    writeReactionInformation(group, context);
}

}
}
}