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
 * @file ReactionCounts.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include "Observable.h"
#include <readdy/common/ParticleTypeTuple.h>
#include <readdy/model/Particle.h>
#include <readdy/model/Context.h>
#include <readdy/model/reactions/Reaction.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class ReactionCounts : public Observable<std::pair<
        std::unordered_map<reactions::Reaction::reaction_id, std::vector<std::size_t>>,
        std::unordered_map<readdy::util::particle_type_pair, std::vector<std::size_t>, readdy::util::particle_type_pair_hasher, readdy::util::particle_type_pair_equal_to>
>> {
public:
    using reaction_counts_order1_map = typename std::tuple_element<0, result_type>::type;
    using reaction_counts_order2_map = typename std::tuple_element<1, result_type>::type;

    ReactionCounts(Kernel* kernel, unsigned int stride);

    ReactionCounts(const ReactionCounts&) = delete;
    ReactionCounts& operator=(const ReactionCounts&) = delete;
    ReactionCounts(ReactionCounts&&) = default;
    ReactionCounts& operator=(ReactionCounts&&) = delete;

    virtual ~ReactionCounts();

    void flush() override;

    /*
     * Initialize the maps corresponding to first and second order reaction counts. If they were not used before, that means creating key-value pairs in
     * the maps and setting the values, which are vectors, to the correct size. If they were used before, all counts within the value-vectors will be
     * filled with zeros. This is used for the reaction-counts object in the state-model as well as the result object of the corresponding observable.
     */
    static void
    initializeCounts(std::pair<ReactionCounts::reaction_counts_order1_map, ReactionCounts::reaction_counts_order2_map> &reactionCounts,
                     const readdy::model::Context &ctx);

protected:
    void initialize(Kernel* kernel) override;

    void initializeDataSet(File &file, const std::string &dataSetName, unsigned int flushStride) override;

    void append() override;

    void assignCountsToResult(const result_type &from, result_type &to);

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
