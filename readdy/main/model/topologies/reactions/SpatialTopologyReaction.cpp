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

#include <readdy/model/topologies/reactions/SpatialTopologyReaction.h>

namespace readdy {
namespace model {
namespace top {
namespace reactions {

SpatialTopologyReaction::SpatialTopologyReaction(const std::string &name, const util::particle_type_pair &types,
                                               const util::particle_type_pair &types_to, const scalar rate,
                                               const scalar radius, STRMode mode)
        : _name(name), _types(types), _types_to(types_to), _rate(rate), _radius(radius), _mode(mode) {}

const std::string &SpatialTopologyReaction::name() const {
    return _name;
}

const particle_type_type SpatialTopologyReaction::type1() const {
    return std::get<0>(_types);
}

const particle_type_type SpatialTopologyReaction::type2() const {
    return std::get<1>(_types);
}

const scalar SpatialTopologyReaction::rate() const {
    return _rate;
}

const scalar SpatialTopologyReaction::radius() const {
    return _radius;
}

const util::particle_type_pair &SpatialTopologyReaction::types() const {
    return _types;
}

const particle_type_type SpatialTopologyReaction::type_to1() const {
    return std::get<0>(_types_to);
}

const particle_type_type SpatialTopologyReaction::type_to2() const {
    return std::get<1>(_types_to);
}

const util::particle_type_pair &SpatialTopologyReaction::types_to() const {
    return _types_to;
}

const bool SpatialTopologyReaction::connect() const {
    return _mode == STRMode::CONNECT || _mode == STRMode::CONNECT_ALLOW_SELF;
}

const bool SpatialTopologyReaction::allow_self_connection() const {
    return _mode == STRMode::CONNECT_ALLOW_SELF;
}


}
}
}
}