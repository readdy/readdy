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
 * Declaration of fusion reactions, i.e., A+B->C. Additionally to the types, they also have a rate, an educt
 * distance and two weights, which determine where (on a line between A and B) C should be placed.
 *
 * @file Fusion.h
 * @brief Header file containing the declaration for fusion reactions, i.e., A+B->C.
 * @author clonker
 * @date 20.06.16
 */

#pragma once
#include "Reaction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

class Fusion : public Reaction {
    using super = Reaction;
public:
    Fusion(const std::string &name, ParticleTypeId from1, ParticleTypeId from2, ParticleTypeId to,
           const scalar rate, const scalar eductDistance, const scalar weight1 = 0.5,
           const scalar weight2 = 0.5) : Reaction(name, rate, eductDistance, 0, 2, 1){
        super::_weight1 = weight1;
        super::_weight2 = weight2;
        _educts = {from1, from2};
        _products = {to};

        const auto sum = weight1 + weight2;
        if (sum != 1) {
            this->_weight1 /= sum;
            this->_weight2 /= sum;
            log::warn("The weights did not add up to 1, they were changed to weight1={}, weight2={}",
                                 this->_weight1, this->_weight2);
        }
    }

    const ParticleTypeId getFrom1() const {
        return _educts[0];
    }

    const ParticleTypeId getFrom2() const {
        return _educts[1];
    }

    const ParticleTypeId getTo() const {
        return _products[0];
    }

    const ReactionType type() const override {
        return ReactionType::Fusion;
    }


};
NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
