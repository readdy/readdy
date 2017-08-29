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

SpatialTopologyReaction::SpatialTopologyReaction(std::string name, util::particle_type_pair types,
                                                 topology_type_pair top_types,
                                                 util::particle_type_pair types_to, topology_type_pair top_types_to,
                                                 scalar rate, scalar radius, STRMode mode)
        : _name(std::move(name)), _types(std::move(types)), _types_to(std::move(types_to)), _rate(rate),
          _radius(radius), _mode(mode), _top_types(std::move(top_types)), _top_types_to(std::move(top_types_to)) {}

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

const bool SpatialTopologyReaction::allow_self_connection() const {
    return _mode == STRMode::TT_FUSION_ALLOW_SELF;
}

const topology_type_type SpatialTopologyReaction::top_type1() const {
    return std::get<0>(_top_types);
}

const topology_type_type SpatialTopologyReaction::top_type2() const {
    return std::get<1>(_top_types);
}

const topology_type_type SpatialTopologyReaction::top_type_to1() const {
    return std::get<0>(_top_types_to);
}

const topology_type_type SpatialTopologyReaction::top_type_to2() const {
    return std::get<1>(_top_types_to);
}

bool SpatialTopologyReaction::is_topology_particle_reaction() const {
    return top_type2() == topology_type_empty;
}

bool SpatialTopologyReaction::is_topology_topology_reaction() const {
    return !is_topology_particle_reaction();
}

const bool SpatialTopologyReaction::is_enzymatic() const {
    return _mode == STRMode::TT_ENZYMATIC || _mode == STRMode::TP_ENZYMATIC;
}

const bool SpatialTopologyReaction::is_fusion() const {
    return _mode == STRMode::TT_FUSION || _mode == STRMode::TT_FUSION_ALLOW_SELF || _mode == STRMode::TP_FUSION;
}

const STRMode &SpatialTopologyReaction::mode() const {
    return _mode;
}

constexpr const char STRParser::arrow[];


}
}
}
}