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
 * Definition and implementation of ParticleTypePair and ParticleTypePair hasher. They serve as a tuple that is order-
 * independent.
 *
 * @file ParticleTypePair.h
 * @brief Header file containing ParticleTypePair and ParticleTypePairHasher.
 * @author clonker
 * @date 03.08.16
 */

#ifndef READDY_MAIN_PARTICLETYPEPAIR_H
#define READDY_MAIN_PARTICLETYPEPAIR_H

#include <cstddef>
#include <tuple>
#include "hash.h"

namespace readdy {
namespace util {

struct ParticleTypePair {
    readdy::model::Particle::type_type t1, t2;

    ParticleTypePair(readdy::model::Particle::type_type t1, readdy::model::Particle::type_type t2) {
        if (t1 <= t2) {
            ParticleTypePair::t1 = t1;
            ParticleTypePair::t2 = t2;
        } else {
            ParticleTypePair::t1 = t2;
            ParticleTypePair::t2 = t1;
        }
    }

    friend std::size_t hash_value(const ParticleTypePair &pair) {
        std::size_t seed = 0;
        hash::combine(seed, pair.t1);
        hash::combine(seed, pair.t2);
        return seed;
    }

    friend bool operator==(const ParticleTypePair &p1, const ParticleTypePair &p2) {
        return p1.t1 == p2.t1 && p1.t2 == p2.t2;
    }
};

class ParticleTypePairHasher {
public:
    std::size_t operator()(const ParticleTypePair &k) const {
        return hash_value(k);
    }

    std::size_t operator()(const std::tuple<readdy::model::Particle::type_type, readdy::model::Particle::type_type> &k) const {
        std::size_t seed = 0;
        const auto &t1 = std::get<0>(k);
        const auto &t2 = std::get<1>(k);
        if (t1 <= t2) {
            hash::combine(seed, t1);
            hash::combine(seed, t2);
        } else {
            hash::combine(seed, t2);
            hash::combine(seed, t1);
        }
        return seed;
    }
};

}
}
#endif //READDY_MAIN_PARTICLETYPEPAIR_H
