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
 * @file ReactionRecord.h
 * @brief << brief description >>
 * @author clonker
 * @date 07.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once
#include <spdlog/fmt/ostr.h>
#include <readdy/common/common.h>
#include <readdy/model/Particle.h>
#include <readdy/model/reactions/Reaction.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

struct ReactionRecord {
    int type {0}; // int corresponding to the enum readdy::model::reactions::ReactionType
    std::array<Particle::id_type, 2> educts {{0, 0}};
    std::array<Particle::id_type, 2> products {{0, 0}};
    std::array<Particle::type_type, 2> types_from {{0, 0}};
    Vec3 where {0, 0, 0};
    // std::size_t reactionIndex {0};
    readdy::model::reactions::Reaction::reaction_id id {0}; // global unique reaction id

    friend std::ostream &operator<<(std::ostream &os, const ReactionRecord &record);
};

using reaction_counts_map = std::unordered_map<reactions::Reaction::reaction_id, std::size_t>;

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
