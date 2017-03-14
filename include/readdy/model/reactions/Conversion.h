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
 * This file contains the declaration for conversion reactions, i.e., A->B. They are assigned to two types and happen
 * with a certain rate.
 *
 * @file Conversion.h
 * @brief Declaration of Conversion reactions, i.e., A->B.
 * @author clonker
 * @date 20.06.16
 */

#pragma once
#include "Reaction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

class Conversion : public Reaction<1> {

public:
    Conversion(const std::string &name, particle_type_type typeFrom, particle_type_type typeTo, const double rate) :
            Reaction(name, rate, 0, 0, 1) {
        educts = {typeFrom};
        products = {typeTo};
    }

    const unsigned int getTypeFrom() const {
        return educts[0];
    }

    const unsigned int getTypeTo() const {
        return products[0];
    }

    virtual const ReactionType getType() override {
        return ReactionType::Conversion;
    }
};
NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
