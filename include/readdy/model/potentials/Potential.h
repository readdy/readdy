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
 * This header file contains the base class for all potentials. In particular, all potentials shall have
 * a uuid by which they can be identified, a name that describes the type and an order, which indicates the number of
 * particles which are interacting with one another.
 *
 * @file Potential.h
 * @brief Header file containing the base class for all potentials.
 * @author clonker
 * @date 31.05.16
 */

#ifndef READDY_MAIN_POTENTIAL_H
#define READDY_MAIN_POTENTIAL_H

#include <string>
#include <array>
#include <readdy/model/Particle.h>

namespace readdy {
namespace model {
namespace potentials {

class Potential {
public:
    using id_t = short;
    
    Potential(const std::string &name, const int order) : name(name), order(order), id(counter++) { }

    virtual ~Potential() = default;

    const id_t getId() const {
        return id;
    }

    const std::string &getName() const {
        return name;
    }

    const int getOrder() const {
        return order;
    }

    virtual double getMaximalForce(double kbt) const noexcept = 0;

private:
    static id_t counter;

    const std::string name;
    const int order;
    const id_t id;
};

}
}
}

#endif //READDY_MAIN_POTENTIAL_H
