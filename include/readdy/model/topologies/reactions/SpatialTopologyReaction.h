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
 * @file TopologyFusionReaction.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/macros.h>
#include <readdy/common/common.h>
#include <readdy/common/ParticleTypeTuple.h>
#include <readdy/model/topologies/TopologyParticleTypeMap.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(reactions)

/**
 * TT - Topology-Topology
 * TP - Topology-Particle
 */
enum class STRMode {
    TT_ENZYMATIC = 0, TT_FUSION, TT_FUSION_ALLOW_SELF, TP_ENZYMATIC, TP_FUSION
};

class SpatialTopologyReaction {
public:
    SpatialTopologyReaction(std::string name, util::particle_type_pair types, topology_type_pair top_types,
                                util::particle_type_pair types_to, topology_type_pair top_types_to, scalar rate,
                                scalar radius, STRMode mode);

    ~SpatialTopologyReaction() = default;

    SpatialTopologyReaction(const SpatialTopologyReaction &) = default;

    SpatialTopologyReaction &operator=(const SpatialTopologyReaction &) = default;

    SpatialTopologyReaction(SpatialTopologyReaction &&) = default;

    SpatialTopologyReaction &operator=(SpatialTopologyReaction &&) = default;

    const std::string &name() const;

    const particle_type_type type1() const;

    const particle_type_type type2() const;

    const util::particle_type_pair &types() const;

    const particle_type_type type_to1() const;

    const particle_type_type type_to2() const;

    const util::particle_type_pair &types_to() const;

    const topology_type_type top_type1() const;

    const topology_type_type top_type2() const;

    const topology_type_type top_type_to1() const;

    const topology_type_type top_type_to2() const;

    bool is_topology_particle_reaction() const;

    bool is_topology_topology_reaction() const;

    const bool is_enzymatic() const;

    const bool is_fusion() const;

    const scalar rate() const;

    const scalar radius() const;

    const bool allow_self_connection() const;

    const STRMode &mode() const;

private:
    std::string _name;
    util::particle_type_pair _types;
    util::particle_type_pair _types_to;
    topology_type_pair _top_types;
    topology_type_pair _top_types_to;
    scalar _rate;
    scalar _radius;
    STRMode _mode;
};

class STRParser {
public:
    /**
     *  Pass descriptor of form for Topology<->Topology
     * - fusion type:
     *     name: T1 (p1) + T2 (p2) -> T3 (p3 + p4) [self=true|false]
     * - enzymatic type:
     *     name: T1 (p1) + T2 (p2) -> T3 (p3) + T4 (p4)
     *
     * and for Topology<->Particle
     * - fusion type:
     *     name: T1 (p1) + p2 -> T2 (p3 + p4)
     * - enzymatic type:
     *     name: T1 (p1) + p2 -> T2 (p3) + p4
     *
     * @param descriptor
     * @param rate
     * @param radius
     * @return
     */
    SpatialTopologyReaction parse(const std::string& descriptor, scalar rate, scalar radius) const {
        log::trace("begin parsing \"{}\"", descriptor);
        auto arrowPos = descriptor.find(arrow);
        if(arrowPos == descriptor.npos) {
            throw std::invalid_argument(fmt::format(
                    "the descriptor must contain an arrow (\"{}\") to indicate lhs and rhs.", arrow
            ));
        }
        if(descriptor.find(arrow, arrowPos+1) != descriptor.npos) {
            throw std::invalid_argument(fmt::format(
                    "the descriptor must not contain more than one arrow (\"{}\").", arrow
            ));
        }
        auto lhs = descriptor.substr(0, arrowPos);
        auto rhs = descriptor.substr(arrowPos + std::strlen(arrow), descriptor.npos);
        log::debug("got lhs={}, rhs={}", lhs, rhs);
    }

private:
    static constexpr const char arrow[] = "->";
};

NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
