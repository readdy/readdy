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

#include <readdy/common/macros.h>
#include <readdy/model/Particle.h>
#include <readdy/common/optional.h>
#include <readdy/model/observables/Observable.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

enum class ReactionType {
    DECAY = 0, CONVERSION, FUSION, FISSION, ENZYMATIC
};

struct ReactionRecord {
    ReactionType type;
    observables::time_step_type when;
    Particle::id_type educts[2];
    Particle::id_type products[2];
    std::experimental::optional<Vec3> where;
};

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
#endif //READDY_MAIN_REACTIONRECORD_H
