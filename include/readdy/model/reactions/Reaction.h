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
 * Base class for all possible types of reaction. Currently:
 *  - Conversion: A -> B
 *  - Enzymatic: A + C -> B + C where C is a catalyst
 *  - Fission: A -> B + C
 *  - Fusion: A + B -> C
 *
 * @file Reactions.h
 * @brief Reaction base class.
 * @author clonker
 * @date 17.06.16
 */

#pragma once
#include <string>
#include <ostream>
#include <utility>
#include <spdlog/fmt/ostr.h>
#include <readdy/model/Particle.h>
#include <readdy/model/RandomProvider.h>
#include <readdy/common/make_unique.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

enum class ReactionType { Conversion, Fusion, Fission, Enzymatic, Decay };

std::ostream& operator<<(std::ostream& os, const ReactionType& reactionType);

class Reaction {
public:
    using rnd_normal = std::function<Vec3(const scalar, const scalar)>;
    using reaction_id = unsigned short;

    Reaction(std::string name, scalar rate, scalar eductDistance, scalar productDistance, std::uint8_t nEducts,
             std::uint8_t nProducts);

    virtual ~Reaction() = default;

    virtual const ReactionType type() const = 0;

    const std::string &name() const;

    const reaction_id id() const;

    const scalar rate() const;

    std::uint8_t nEducts() const;

    std::uint8_t nProducts() const;

    const scalar eductDistance() const;

    const scalar eductDistanceSquared() const;

    const scalar productDistance() const;

    friend std::ostream &operator<<(std::ostream &os, const Reaction &reaction);

    const std::array<particle_type_type, 2> &educts() const;

    const std::array<particle_type_type, 2> &products() const;

    const scalar weight1() const;

    const scalar weight2() const;


protected:
    static reaction_id counter;

    std::uint8_t _nEducts;
    std::uint8_t _nProducts;
    std::array<particle_type_type, 2> _educts {{0, 0}};
    std::array<particle_type_type, 2> _products {{0, 0}};
    std::string _name;
    reaction_id _id;
    scalar _rate;
    scalar _eductDistance, _eductDistanceSquared;
    scalar _productDistance;

    scalar _weight1 = .5, _weight2 = .5;
};

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
