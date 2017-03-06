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

#ifndef READDY_MAIN_REACTION_H
#define READDY_MAIN_REACTION_H

#include <string>
#include <ostream>
#include <spdlog/fmt/ostr.h>
#include <readdy/model/Particle.h>
#include <readdy/model/RandomProvider.h>
#include <readdy/common/make_unique.h>

namespace readdy {
namespace model {
namespace reactions {

template<unsigned int N_EDUCTS>
class Reaction {
protected:
    static short counter;
    using particle_type_type = readdy::model::Particle::type_type;
public:
    enum ReactionType { Conversion, Fusion, Fission, Enzymatic, Decay };

    using rnd_normal = std::function<Vec3(const double, const double)>;
    // static constexpr unsigned int n_educts = N_EDUCTS;

    Reaction(const std::string &name, const double rate, const double eductDistance,
             const double productDistance, const unsigned int n_products) :
            name(name),
            id(counter++),
            rate(rate),
            eductDistance(eductDistance),
            eductDistanceSquared(eductDistance * eductDistance),
            productDistance(productDistance),
            _n_products(n_products) {}

    virtual ~Reaction() = default;

    virtual const ReactionType getType() = 0;

    const std::string &getName() const {
        return name;
    }

    const short getId() const {
        return id;
    }

    const double getRate() const {
        return rate;
    }

    const unsigned int getNEducts() const {
        return _n_educts;
    }

    const unsigned int getNProducts() const {
        return _n_products;
    }

    const double getEductDistance() const {
        return eductDistance;
    }

    const double getEductDistanceSquared() const {
        return eductDistanceSquared;
    }

    const double getProductDistance() const {
        return productDistance;
    }

    friend std::ostream &operator<<(std::ostream &os, const Reaction &reaction) {
        os << "Reaction(\"" << reaction.name << "\", N_Educts=" << reaction._n_educts << ", N_Products="
           << reaction._n_products << ", (";
        for (unsigned int i = 0; i < reaction._n_educts; i++) {
            if (i > 0) os << ",";
            os << reaction.educts[i];
        }
        os << ") -> (";
        for (unsigned int i = 0; i < reaction._n_products; i++) {
            if (i > 0) os << ",";
            os << reaction.products[i];
        }
        os << "), rate=" << reaction.rate << ", eductDist=" << reaction.eductDistance << ", prodDist="
           << reaction.productDistance << ")";
        return os;
    }

    const std::array<particle_type_type, N_EDUCTS> &getEducts() const {
        return educts;
    }

    const std::array<particle_type_type, 2> &getProducts() const {
        return products;
    }

    Reaction(const Reaction &rhs)
            : _n_educts(rhs._n_educts), _n_products(rhs._n_products), educts(rhs.educts),
              products(rhs.products), name(rhs.name), id(rhs.id), rate(rhs.rate),
              eductDistance(rhs.eductDistance), productDistance(rhs.productDistance),
              eductDistanceSquared(rhs.eductDistanceSquared) {
    }

    const double getWeight1() const {
        return weight1;
    }

    const double getWeight2() const {
        return weight2;
    }


protected:
    const unsigned int _n_educts = N_EDUCTS;
    const unsigned int _n_products;
    std::array<particle_type_type, N_EDUCTS> educts;
    std::array<particle_type_type, 2> products;
    const std::string name;
    const short id;
    const double rate;
    const double eductDistance, eductDistanceSquared;
    const double productDistance;

    double weight1 = .5, weight2 = .5;
};

}
}
}

#endif //READDY_MAIN_REACTION_H
