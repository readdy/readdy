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
 * @file Reactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.10.16
 */

#include "readdy/model/reactions/Reaction.h"
#include "readdy/model/reactions/ReactionRecord.h"

namespace readdy {
namespace model {
namespace reactions {

std::ostream& operator<<(std::ostream& os, const ReactionType& reactionType) {
    switch (reactionType) {
        case ReactionType::Decay: os << "Decay"; break;
        case ReactionType::Conversion: os << "Conversion"; break;
        case ReactionType::Fusion: os << "Fusion"; break;
        case ReactionType::Fission: os << "Fission"; break;
        case ReactionType::Enzymatic: os << "Enzymatic"; break;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const ReactionRecord& record) {
    auto type = ReactionType(record.type);
    os << "ReactionRecord[type: " << type;
    switch (type) {
        case ReactionType::Decay:{
            os << ", educt: " << record.educts[0];
            break;
        }
        case ReactionType::Conversion: {
            os << ", educt: " << record.educts[0] << ", product: " << record.products[0];
            break;
        }
        case ReactionType::Fusion: {
            os << ", educts: [" << record.educts[0] << "," << record.educts[1] << "], product: " << record.products[0];
            break;
        }
        case ReactionType::Fission: {
            os << ", educt: " << record.educts[0] << ", products: [" << record.products[0] << "," << record.products[1] << "]";
            break;
        }
        case ReactionType::Enzymatic: {
            os << ", educts: [" << record.educts[0] << "," << record.educts[1] << "]";
            os << ", products: [" << record.products[0] << "," << record.products[1] << "]";
            break;
        }
    }
    os << ", location: " << record.where << "]";
    return os;
}

Reaction::Reaction(std::string name, const scalar rate, const scalar eductDistance, const scalar productDistance,
                   std::uint8_t nEducts, std::uint8_t nProducts)
        : _name(std::move(name)), _id(counter++), _rate(rate), _eductDistance(eductDistance),
          _eductDistanceSquared(eductDistance * eductDistance), _productDistance(productDistance),
          _nEducts(nEducts), _nProducts(nProducts) {}

const std::string &Reaction::name() const {
    return _name;
}

const Reaction::reaction_id Reaction::id() const {
    return _id;
}

const scalar Reaction::rate() const {
    return _rate;
}

std::uint8_t Reaction::nEducts() const {
    return _nEducts;
}

std::uint8_t Reaction::nProducts() const {
    return _nProducts;
}

const scalar Reaction::eductDistance() const {
    return _eductDistance;
}

const scalar Reaction::eductDistanceSquared() const {
    return _eductDistanceSquared;
}

const scalar Reaction::productDistance() const {
    return _productDistance;
}

std::ostream &operator<<(std::ostream &os, const Reaction &reaction) {
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

const std::array<particle_type_type, 2> &Reaction::educts() const {
    return _educts;
}

const std::array<particle_type_type, 2> &Reaction::products() const {
    return _products;
}

const scalar Reaction::weight1() const {
    return _weight1;
}

const scalar Reaction::weight2() const {
    return _weight2;
}

Reaction::reaction_id Reaction::counter = 0;

}
}
}
