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
 * Decay reactions simply remove particles of a certain type with a certain rate of the system.
 *
 * @file Death.h
 * @brief This file contains the declaration of decay reactions.
 * @author clonker
 * @date 21.06.16
 */


#pragma once
#include "Reaction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

class Decay : public Reaction {

public:
    Decay(const std::string &name, ParticleTypeId typeFrom, const scalar rate) : Reaction(name, rate, 0, 0, 1, 0) {
        _educts[0] = typeFrom;
    }

    const ParticleTypeId getTypeFrom() const { return _educts[0]; }

    const ReactionType type() const override { return ReactionType::Decay; }
};
NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
