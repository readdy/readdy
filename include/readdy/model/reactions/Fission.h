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

#ifndef READDY_MAIN_FISSION_H
#define READDY_MAIN_FISSION_H

#include <readdy/common/logging.h>
#include "Reaction.h"

namespace readdy {
namespace model {
namespace reactions {

class Fission : public Reaction<1> {
    using super = Reaction<1>;
public:
    Fission(const std::string &name, particle_type_type from, particle_type_type to1, particle_type_type to2,
            const double rate, const double productDistance, const double weight1 = 0.5,
            const double weight2 = 0.5) :
            Reaction(name, rate, 0, productDistance, 2) {
        super::weight1 = weight1;
        super::weight2 = weight2;
        educts = {from};
        products = {to1, to2};
        const auto sum = weight1 + weight2;
        if (sum != 1) {
            this->weight1 /= sum;
            this->weight2 /= sum;
            log::warn("The weights did not add up to 1, they were changed to weight1={}, weight2={}",
                                         this->weight1, this->weight2);
        }
    }

    const unsigned int getFrom() const {
        return educts[0];
    }

    const unsigned int getTo1() const {
        return products[0];
    }

    const unsigned int getTo2() const {
        return products[1];
    }

    virtual const ReactionType getType() override {
        return ReactionType::Fission;
    }

};
}
}
}
#endif //READDY_MAIN_FISSION_H
