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
 * Write datasets that describe the context of the system that was used.
 * Since these datasets are meta-data, turn off compression, i.e. no blosc-filter is necessary for reading this.
 *
 * @file IOUtils.cpp
 * @brief Utils for writing to files that concern the model, e.g. writing reaction information.
 * @author clonker
 * @author chrisfroe
 * @date 10.03.17
 * @copyright GPL-3
 */

#include <json.hpp>
#include <readdy/model/IOUtils.h>

using json = nlohmann::json;

namespace readdy {
namespace model {
namespace ioutils {

void writeReactionInformation(h5rd::Group &group, const Context &context) {
    const auto &reactionRegistry = context.reactions();
    std::vector<ReactionInfo> reactionInfos;
    for (const auto &reaction : reactionRegistry.order1Flat()) {
        const auto &r = reaction;
        ReactionInfo info{r->name().c_str(), r->id(), r->nEducts(), r->nProducts(),
                          r->rate(), r->eductDistance(),
                          r->productDistance(), r->educts(), r->products()};
        reactionInfos.push_back(info);
    }
    for (const auto &reaction : reactionRegistry.order2Flat()) {
        const auto &r = reaction;
        ReactionInfo info{r->name().c_str(), r->id(), r->nEducts(), r->nProducts(),
                          r->rate(), r->eductDistance(),
                          r->productDistance(), r->educts(), r->products()};
        reactionInfos.push_back(info);
    }
    if(!reactionInfos.empty()) {
        auto types = getReactionInfoMemoryType(group.parentFile());
        h5rd::dimensions dims = {h5rd::UNLIMITED_DIMS};
        h5rd::dimensions extent = {reactionInfos.size()};
        std::vector<h5rd::Filter*> filters;
        auto infoDataSet = group.createDataSet("registered_reactions", extent, dims, std::get<0>(types),
                                               std::get<1>(types), filters);
        infoDataSet->append(extent, reactionInfos.data());
    }
}

void writeGeneralContextInformation(h5rd::Group &group, const Context &context) {
    json j;
    j["kbt"] = context.kBT();
    j["box_volume"] = context.boxVolume();
    j["box_size"] = context.boxSize();
    j["pbc"] = context.periodicBoundaryConditions();
    group.write("general", j.dump());
}

void writeParticleTypeInformation(h5rd::Group &group, const Context &context) {
    auto h5types = getParticleTypeInfoType(group.parentFile());

    const auto &types = context.particleTypes().typeMapping();
    std::vector<ParticleTypeInfo> type_info_vec;
    for (const auto &p_type : types) {
        ParticleTypeInfo info{p_type.first.c_str(), p_type.second,
                              context.particleTypes().diffusionConstantOf(p_type.first)};
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
            .insert<decltype(std::declval<ReactionInfo>().id)>("id", offsetof(ReactionInfo, id))
            .insert<decltype(std::declval<ReactionInfo>().n_educts)>("n_educts", offsetof(ReactionInfo, n_educts))
            .insert<decltype(std::declval<ReactionInfo>().n_products)>("n_products", offsetof(ReactionInfo, n_products))
            .insert<decltype(std::declval<ReactionInfo>().rate)>("rate", offsetof(ReactionInfo, rate))
            .insert<decltype(std::declval<ReactionInfo>().educt_distance)>("educt_distance", offsetof(ReactionInfo, educt_distance))
            .insert<decltype(std::declval<ReactionInfo>().product_distance)>("product_distance", offsetof(ReactionInfo, product_distance))
            .insertArray<ParticleTypeId, 2>("educt_types", offsetof(ReactionInfo, educt_types))
            .insertArray<ParticleTypeId, 2>("product_types", offsetof(ReactionInfo, product_types))
            .build();
    STDCompoundType sct (nct);
    return std::make_tuple(nct, sct);
};

void writeSimulationSetup(h5rd::Group &group, const Context &context) {
    writeParticleTypeInformation(group, context);
    writeReactionInformation(group, context);
    writeGeneralContextInformation(group, context);
}

}
}
}