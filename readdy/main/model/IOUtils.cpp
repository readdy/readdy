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

namespace readdy {
namespace model {
namespace ioutils {

// @todo save the registered_reactions as a flat vector of ReactionInfo
// @todo omit the reaction index in ReactionInfo, because it is an implementation detail
void writeReactionInformation(h5rd::Group &group, const Context &context) {
    auto subgroup = group.createGroup("./registered_reactions");
    // order1
    const auto &order1_reactions = context.reactions().order1Flat();
    auto n_reactions = order1_reactions.size();
    auto types = getReactionInfoMemoryType(group.parentFile());
    if (n_reactions > 0) {
        std::vector<ReactionInfo> order1_info;
        for (const auto &r : order1_reactions) {
            const auto &reactions_current_type = context.reactions().order1ByType(r->educts()[0]);
            auto it = std::find_if(reactions_current_type.begin(), reactions_current_type.end(),
                                   [&r](const reactions::Reaction *x) { return x->id() == r->id(); });
            if (it != reactions_current_type.end()) {
                std::size_t index = static_cast<std::size_t>(it - reactions_current_type.begin());
                const std::array<particle_type_type, 2> educts = {r->educts()[0], 0};
                ReactionInfo info{r->name().c_str(), index, r->id(), r->nEducts(), r->nProducts(),
                                  r->rate(), r->eductDistance(),
                                  r->productDistance(), educts, r->products()};
                order1_info.push_back(info);
            }
        }
        h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
        h5rd::dimensions extent = {n_reactions};
        std::vector<h5rd::Filter*> filters;
        auto order1_reaction_dset = subgroup.createDataSet("order1_reactions", extent, dims, std::get<0>(types),
                                                           std::get<1>(types), filters);
        order1_reaction_dset->append(extent, order1_info.data());
    }
    // order2
    const auto &order2_reactions = context.reactions().order2Flat();
    n_reactions = order2_reactions.size();
    if (n_reactions > 0) {
        std::vector<ReactionInfo> order2_info;
        for (const auto &r : order2_reactions) {
            const auto &reactions_current_type = context.reactions().order2ByType(r->educts()[0],
                                                                                  r->educts()[1]);
            auto it = std::find_if(reactions_current_type.begin(), reactions_current_type.end(),
                                   [&r](const reactions::Reaction *x) { return x->id() == r->id(); });
            if (it != reactions_current_type.end()) {
                std::size_t index = static_cast<std::size_t>(it - reactions_current_type.begin());
                ReactionInfo info{r->name().c_str(), index, r->id(), r->nEducts(), r->nProducts(),
                                  r->rate(), r->eductDistance(),
                                  r->productDistance(), r->educts(), r->products()};
                order2_info.push_back(info);
            }
        }
        h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
        h5rd::dimensions extent = {n_reactions};
        std::vector<h5rd::Filter*> filters;
        auto order2_reaction_dset = subgroup.createDataSet("order2_reactions", extent, dims, std::get<0>(types),
                                                           std::get<1>(types), filters);
        order2_reaction_dset->append(extent, order2_info.data());
    }
}

void writeParticleTypeInformation(h5rd::Group &group, const Context &context) {
    auto h5types = getParticleTypeInfoType(group.parentFile());

    const auto &types = context.particle_types().typeMapping();
    std::vector<ParticleTypeInfo> type_info_vec;
    for (const auto &p_type : types) {
        ParticleTypeInfo info{p_type.first.c_str(), p_type.second,
                              context.particle_types().diffusionConstantOf(p_type.first)};
        type_info_vec.push_back(info);
    }
    if (!type_info_vec.empty()) {
        h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
        h5rd::dimensions extent = {type_info_vec.size()};
        std::vector<h5rd::Filter*> filters;
        auto dset = group.createDataSet("particle_types", extent, dims, std::get<0>(h5types),
                                        std::get<1>(h5types), filters);
        dset->append(extent, type_info_vec.data());
    }
}

std::tuple<h5rd::NativeCompoundType, h5rd::STDCompoundType> getParticleTypeInfoType(h5rd::Object::ParentFileRef ref) {
    using namespace h5rd;
    NativeCompoundType nct = NativeCompoundTypeBuilder(sizeof(ParticleTypeInfo), std::move(ref))
            .insertString("name", offsetof(ParticleTypeInfo, name))
            .insert<decltype(std::declval<ParticleTypeInfo>().type_id)>("type_id", offsetof(ParticleTypeInfo, type_id))
            .insert<decltype(std::declval<ParticleTypeInfo>().diffusion_constant)>("diffusion_constant", offsetof(ParticleTypeInfo, diffusion_constant))
            .build();
    return std::make_tuple(nct, STDCompoundType(nct));
};

std::tuple<h5rd::NativeCompoundType, h5rd::STDCompoundType> getReactionInfoMemoryType(h5rd::Object::ParentFileRef ref) {
    using namespace h5rd;
    NativeCompoundType nct  = NativeCompoundTypeBuilder(sizeof(ReactionInfo), std::move(ref))
            .insertString("name", offsetof(ReactionInfo, name))
            .insert<decltype(std::declval<ReactionInfo>().index)>("index", offsetof(ReactionInfo, index))
            .insert<decltype(std::declval<ReactionInfo>().id)>("id", offsetof(ReactionInfo, id))
            .insert<decltype(std::declval<ReactionInfo>().n_educts)>("n_educts", offsetof(ReactionInfo, n_educts))
            .insert<decltype(std::declval<ReactionInfo>().n_products)>("n_products", offsetof(ReactionInfo, n_products))
            .insert<decltype(std::declval<ReactionInfo>().rate)>("rate", offsetof(ReactionInfo, rate))
            .insert<decltype(std::declval<ReactionInfo>().educt_distance)>("educt_distance", offsetof(ReactionInfo, educt_distance))
            .insert<decltype(std::declval<ReactionInfo>().product_distance)>("product_distance", offsetof(ReactionInfo, product_distance))
            .insertArray<particle_type_type, 2>("educt_types", offsetof(ReactionInfo, educt_types))
            .insertArray<particle_type_type, 2>("product_types", offsetof(ReactionInfo, product_types))
            .build();
    STDCompoundType sct (nct);
    return std::make_tuple(nct, sct);
};

void writeSimulationSetup(h5rd::Group &group, const Context &context) {
    writeParticleTypeInformation(group, context);
    writeReactionInformation(group, context);
}

}
}
}