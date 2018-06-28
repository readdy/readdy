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
 * This file contains the declaration for enzymatic reactions, i.e., A+C->B+C.
 *
 * @file Enzymatic.h
 * @brief Declaration file of enzymatic reactions, i.e., A+C->B+C.
 * @author clonker
 * @date 20.06.16
 */

#pragma once
#include "Reaction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

class Enzymatic : public Reaction {

public:
    Enzymatic(const std::string &name, ParticleTypeId catalyst, ParticleTypeId from, ParticleTypeId to,
              const scalar rate, const scalar eductDistance) :
            Reaction(name, rate, eductDistance, 0, 2, 2) {
        _educts = {from, catalyst};
        _products = {to, catalyst};
    }


    const ParticleTypeId getCatalyst() const { return _educts[1]; }

    const ParticleTypeId getFrom() const { return _educts[0]; }

    const ParticleTypeId getTo() const { return _products[0]; }

    const ReactionType type() const override { return ReactionType::Enzymatic; }
};
NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
