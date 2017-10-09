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
 * This header file contains the declaration of fission reactions, i.e., A->B+C. The products can be placed at a
 * specified distance and weighted by two members.
 *
 * @file Fission.h
 * @brief Declaration of fission reactions, i.e., A->B+C.
 * @author clonker
 * @date 20.06.16
 */

#pragma once

#include <readdy/common/logging.h>
#include "Reaction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

class Fission : public Reaction {
    using super = Reaction;
public:
    Fission(const std::string &name, particle_type_type from, particle_type_type to1, particle_type_type to2,
            const scalar rate, const scalar productDistance, const scalar weight1 = 0.5,
            const scalar weight2 = 0.5) :
            Reaction(name, rate, 0, productDistance, 1, 2) {
        super::_weight1 = weight1;
        super::_weight2 = weight2;
        _educts = {from};
        _products = {to1, to2};
        const auto sum = weight1 + weight2;
        if (sum != 1) {
            this->_weight1 /= sum;
            this->_weight2 /= sum;
            log::warn("The weights did not add up to 1, they were changed to weight1={}, weight2={}",
                      this->_weight1, this->_weight2);
        }
    }

    const particle_type_type getFrom() const {
        return _educts[0];
    }

    const particle_type_type getTo1() const {
        return _products[0];
    }

    const particle_type_type getTo2() const {
        return _products[1];
    }

    virtual const ReactionType type() const override {
        return ReactionType::Fission;
    }

};
NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
