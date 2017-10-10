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

#pragma once

#include <string>
#include <array>
#include <readdy/model/Particle.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(potentials)

class Potential {
public:
    using id_type = short;
    using order_type = std::uint8_t;

    explicit Potential(order_type order) : _order(order), _id(counter++) { }

    Potential(const Potential&) = default;
    Potential& operator=(const Potential&) = delete;
    Potential(Potential&&) = default;
    Potential& operator=(Potential&& ) = default;

    virtual ~Potential() = default;

    const id_type getId() const {
        return _id;
    }

    order_type getOrder() const {
        return _order;
    }

    virtual scalar getMaximalForce(scalar kbt) const noexcept = 0;

    virtual std::string describe() const = 0;

    virtual std::string type() const = 0;

    friend std::ostream &operator<<(std::ostream &os, const Potential &potential) {
        os << potential.describe();
        return os;
    }

private:
    static id_type counter;

    order_type _order;
    id_type _id;
};

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
