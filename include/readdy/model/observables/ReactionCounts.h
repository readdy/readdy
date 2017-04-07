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
#include <readdy/model/KernelContext.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class ReactionCounts : public Observable<std::pair<
        std::unordered_map<particle_type_type, std::vector<std::size_t>>,
        std::unordered_map<readdy::util::particle_type_pair, std::vector<std::size_t>, readdy::util::particle_type_pair_hasher, readdy::util::particle_type_pair_equal_to>
>> {
public:
    using reaction_counts_order1_map = typename std::tuple_element<0, result_t>::type;
    using reaction_counts_order2_map = typename std::tuple_element<1, result_t>::type;

    ReactionCounts(Kernel *const kernel, unsigned int stride);

    virtual ~ReactionCounts();

    virtual void flush() override;

    /*
     * Initialize the maps corresponding to first and second order reaction counts. If they were not used before, that means creating key-value pairs in
     * the maps and setting the values, which are vectors, to the correct size. If they were used before, all counts within the value-vectors will be
     * filled with zeros. This is used for the reaction-counts object in the state-model as well as the result object of the corresponding observable.
     */
    static void
    initializeCounts(std::pair<ReactionCounts::reaction_counts_order1_map, ReactionCounts::reaction_counts_order2_map> &reactionCounts,
                     const readdy::model::KernelContext &ctx) {
        auto &order1Counts = std::get<0>(reactionCounts);
        auto &order2Counts = std::get<1>(reactionCounts);
        for (const auto &entry1 : ctx.particle_types().type_mapping()) {
            const auto &pType1 = entry1.second;
            const auto numberReactionsOrder1 = ctx.reactions().order1_by_type(pType1).size();
            if (numberReactionsOrder1 > 0) {
                // will create an entry for pType1 if necessary
                auto &countsForType = order1Counts[pType1];
                if (countsForType.empty()) {
                    countsForType.resize(numberReactionsOrder1);
                } else {
                    std::fill(countsForType.begin(), countsForType.end(), 0);
                }
            }
            for (const auto &entry2: ctx.particle_types().type_mapping()) {
                const auto &pType2 = entry2.second;
                if (pType2 < pType1) continue;
                const auto numberReactionsOrder2 = ctx.reactions().order2_by_type(pType1, pType2).size();
                if (numberReactionsOrder2 > 0) {
                    // will create an entry for particle-type-pair if necessary
                    auto &countsForPair = order2Counts[std::tie(pType1, pType2)];
                    if (countsForPair.empty()) {
                        countsForPair.resize(numberReactionsOrder2);
                    } else {
                        std::fill(countsForPair.begin(), countsForPair.end(), 0);
                    }
                }
            }
        }
    }

protected:
    virtual void initialize(Kernel *const kernel) override;

    virtual void initializeDataSet(io::File &file, const std::string &dataSetName, unsigned int flushStride) override;

    virtual void append() override;

    template<typename counts_map_t, typename dset_map_t>
    void writeCountsToDataSets(const counts_map_t &countsMap, dset_map_t &dsetMap) {
        for (const auto &entry : countsMap) {
            auto &dataSet = dsetMap.at(entry.first);
            const auto &counts = entry.second;
            dataSet->append({1, counts.size()}, counts.data());
        }
    }

    template<typename counts_map_t>
    void assignVectorsOfMap(const counts_map_t &from, counts_map_t &to) {
        for (const auto &entry: from) {
            const auto &fromVector = entry.second;
            auto &toVector = to.at(entry.first);
            toVector.assign(fromVector.begin(), fromVector.end());
        }
    }

    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
