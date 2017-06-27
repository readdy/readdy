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
 * @file TopologyFusionReaction.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/reactions/TopologyFusionReaction.h>

namespace readdy {
namespace model {
namespace top {
namespace reactions {

TopologyFusionReaction::TopologyFusionReaction(const std::string &name, const util::particle_type_pair &types,
                                               const util::particle_type_pair &types_to, const scalar rate,
                                               const scalar radius)
        : _name(name), _types(types), _types_to(types_to), _rate(rate), _radius(radius) {}

const std::string &TopologyFusionReaction::name() const {
    return _name;
}

const particle_type_type TopologyFusionReaction::type1() const {
    return std::get<0>(_types);
}

const particle_type_type TopologyFusionReaction::type2() const {
    return std::get<1>(_types);
}

const scalar TopologyFusionReaction::rate() const {
    return _rate;
}

const scalar TopologyFusionReaction::radius() const {
    return _radius;
}

const util::particle_type_pair &TopologyFusionReaction::types() const {
    return _types;
}

const particle_type_type TopologyFusionReaction::type_to1() const {
    return std::get<0>(_types_to);
}

const particle_type_type TopologyFusionReaction::type_to2() const {
    return std::get<1>(_types_to);
}

const util::particle_type_pair &TopologyFusionReaction::types_to() const {
    return _types_to;
}


}
}
}
}