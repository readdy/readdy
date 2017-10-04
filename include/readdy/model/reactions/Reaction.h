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

template<unsigned int N_EDUCTS>
class Reaction {
protected:
    static short counter;
public:
    using rnd_normal = std::function<Vec3(const scalar, const scalar)>;

    Reaction(std::string name, const scalar rate, const scalar eductDistance,
             const scalar productDistance, const unsigned int nProducts) :
            _name(std::move(name)),
            _id(counter++),
            _rate(rate),
            _eductDistance(eductDistance),
            _eductDistanceSquared(eductDistance * eductDistance),
            _productDistance(productDistance),
            _nProducts(nProducts) {}

    virtual ~Reaction() = default;

    virtual const ReactionType type() = 0;

    const std::string &name() const {
        return _name;
    }

    const short id() const {
        return _id;
    }

    const scalar rate() const {
        return _rate;
    }

    const unsigned int nEducts() const {
        return _nEducts;
    }

    const unsigned int nProducts() const {
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
        os << "Reaction(\"" << reaction._name << "\", N_Educts=" << reaction._nEducts << ", N_Products="
           << reaction._nProducts << ", (";
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

    const std::array<particle_type_type, N_EDUCTS> &educts() const {
        return _educts;
    }

    const std::array<particle_type_type, 2> &products() const {
        return _products;
    }

    const scalar weight1() const {
        return _weight1;
    }

    const scalar weight2() const {
        return _weight2;
    }


protected:
    unsigned int _nEducts = N_EDUCTS;
    unsigned int _nProducts;
    std::array<particle_type_type, N_EDUCTS> _educts;
    std::array<particle_type_type, 2> _products {{0, 0}};
    std::string _name;
    short _id;
    scalar _rate;
    scalar _eductDistance, _eductDistanceSquared;
    scalar _productDistance;

    scalar _weight1 = .5, _weight2 = .5;
};

template<unsigned int N> short Reaction<N>::counter = 0;

NAMESPACE_END(reactions)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
