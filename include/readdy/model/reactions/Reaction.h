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
#include <memory>
#include <spdlog/fmt/ostr.h>
#include <readdy/model/Particle.h>
#include <readdy/model/RandomProvider.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(reactions)

enum class ReactionType { Conversion, Fusion, Fission, Enzymatic, Decay };

inline std::ostream& operator<<(std::ostream& os, const ReactionType& reactionType) {
    switch (reactionType) {
        case ReactionType::Decay: os << "Decay"; break;
        case ReactionType::Conversion: os << "Conversion"; break;
        case ReactionType::Fusion: os << "Fusion"; break;
        case ReactionType::Fission: os << "Fission"; break;
        case ReactionType::Enzymatic: os << "Enzymatic"; break;
    }
    return os;
}

class Reaction {
public:
    using ReactionId = unsigned short;

    Reaction(std::string name, scalar rate, scalar eductDistance, scalar productDistance, std::uint8_t nEducts,
             std::uint8_t nProducts)
            : _name(std::move(name)), _id(counter++), _rate(rate), _eductDistance(eductDistance),
              _eductDistanceSquared(eductDistance * eductDistance), _productDistance(productDistance),
              _nEducts(nEducts), _nProducts(nProducts) {};

    virtual ~Reaction() = default;

    virtual const ReactionType type() const = 0;

    const std::string &name() const {
        return _name;
    }

    const ReactionId id() const {
        return _id;
    }

    const scalar rate() const {
        return _rate;
    }

    std::uint8_t nEducts() const {
        return _nEducts;
    }

    std::uint8_t nProducts() const {
        return _nProducts;
    }

    const scalar eductDistance() const {
        return _eductDistance;
    }

    const scalar eductDistanceSquared() const {
        return _eductDistanceSquared;
    }

    const scalar productDistance() const {
        return _productDistance;
    }

    friend std::ostream &operator<<(std::ostream &os, const Reaction &reaction) {
        os << "Reaction(\"" << reaction._name << "\", N_Educts=" << std::to_string(reaction._nEducts) << ", N_Products="
           << std::to_string(reaction._nProducts) << ", (";
        for (unsigned int i = 0; i < reaction._nEducts; i++) {
            if (i > 0) os << ",";
            os << reaction._educts[i];
        }
        os << ") -> (";
        for (unsigned int i = 0; i < reaction._nProducts; i++) {
            if (i > 0) os << ",";
            os << reaction._products[i];
        }
        os << "), rate=" << reaction._rate << ", eductDist=" << reaction._eductDistance << ", prodDist="
           << reaction._productDistance << ")";
        return os;
    }

    const std::array<ParticleTypeId, 2> &educts() const {
        return _educts;
    };

    const std::array<ParticleTypeId, 2> &products() const {
        return _products;
    };

    const scalar weight1() const {
        return _weight1;
    }

    const scalar weight2() const {
        return _weight2;
    }


protected:
    static ReactionId counter;

    std::uint8_t _nEducts;
    std::uint8_t _nProducts;
    std::array<ParticleTypeId, 2> _educts {{0, 0}};
    std::array<ParticleTypeId, 2> _products {{0, 0}};
    std::string _name;
    ReactionId _id;
    scalar _rate;
    scalar _eductDistance, _eductDistanceSquared;
    scalar _productDistance;

    scalar _weight1 = .5, _weight2 = .5;
};

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
