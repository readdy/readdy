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

#ifndef READDY_MAIN_REACTIONRECORD_H
#define READDY_MAIN_REACTIONRECORD_H

#include <readdy/common/common.h>
#include <readdy/model/Particle.h>
#include <spdlog/fmt/ostr.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

enum class ReactionType {
    DECAY = 0, CONVERSION, FUSION, FISSION, ENZYMATIC
};

inline std::ostream& operator<<(std::ostream& os, const ReactionType& reactionType);

struct ReactionRecord {
    ReactionType type = ReactionType::DECAY;
    time_step_type when = 0;
    Particle::id_type educts[2] = {0, 0};
    Particle::id_type products[2] = {0, 0};
    // 1st element tells if 2nd argument is initialized
    std::tuple<bool, Vec3> where = std::make_tuple(false, Vec3(0,0,0));

    friend std::ostream &operator<<(std::ostream &os, const ReactionRecord &record);
};

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
#endif //READDY_MAIN_REACTIONRECORD_H
