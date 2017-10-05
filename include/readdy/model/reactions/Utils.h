/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          *
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
 * @file Utils.h
 * @brief << brief description >>
 * @author chrisfroe
 * @date 05.10.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/common.h>
#include <readdy/model/Particle.h>
#include <readdy/model/reactions/Reaction.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)
NAMESPACE_BEGIN(utils)

using reaction_counts_map = std::unordered_map<reactions::Reaction::reaction_id, std::size_t>;

inline void zeroReactionCounts(reaction_counts_map &reactionCounts,
                               const std::unordered_map<Reaction::reaction_id, Reaction *> &allReactions) {
    reactionCounts.clear();
    for (const auto &reaction : allReactions) {
        const auto &r = reaction.second;
        reactionCounts[r->id()] = 0;
    }
}

NAMESPACE_END(utils)
NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
