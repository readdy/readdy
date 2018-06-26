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
#include <readdy/common/string.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
class TopologyRegistry;
NAMESPACE_BEGIN(reactions)

/**
 * TT - Topology-Topology
 * TP - Topology-Particle
 */
enum class STRMode {
    TT_ENZYMATIC = 0, TT_FUSION, TT_FUSION_ALLOW_SELF, TP_ENZYMATIC, TP_FUSION
};

class STRParser;

class SpatialTopologyReaction {
public:
    SpatialTopologyReaction(std::string name, util::particle_type_pair types, topology_type_pair top_types,
                            util::particle_type_pair types_to, topology_type_pair top_types_to, scalar rate,
                            scalar radius, STRMode mode)
            : _name(std::move(name)), _types(std::move(types)), _types_to(std::move(types_to)), _rate(rate),
              _radius(radius), _mode(mode), _top_types(std::move(top_types)), _top_types_to(std::move(top_types_to)) {};

    ~SpatialTopologyReaction() = default;

    SpatialTopologyReaction(const SpatialTopologyReaction &) = default;

    SpatialTopologyReaction &operator=(const SpatialTopologyReaction &) = default;

    SpatialTopologyReaction(SpatialTopologyReaction &&) = default;

    SpatialTopologyReaction &operator=(SpatialTopologyReaction &&) = default;

    const std::string &name() const {
        return _name;
    }

    const ParticleTypeId type1() const {
        return std::get<0>(_types);
    }

    const ParticleTypeId type2() const {
        return std::get<1>(_types);
    }

    const util::particle_type_pair &types() const {
        return _types;
    }

    const ParticleTypeId type_to1() const {
        return std::get<0>(_types_to);
    }

    const ParticleTypeId type_to2() const {
        return std::get<1>(_types_to);
    }

    const util::particle_type_pair &types_to() const {
        return _types_to;
    }

    const TopologyTypeId top_type1() const {
        return std::get<0>(_top_types);
    }

    const TopologyTypeId top_type2() const {
        return std::get<1>(_top_types);
    }

    const TopologyTypeId top_type_to1() const {
        return std::get<0>(_top_types_to);
    }

    const TopologyTypeId top_type_to2() const {
        return std::get<1>(_top_types_to);
    }

    bool is_topology_particle_reaction() const {
        return top_type2() == EmptyTopologyId;
    }

    bool is_topology_topology_reaction() const {
        return !is_topology_particle_reaction();
    }

    const bool is_enzymatic() const {
        return _mode == STRMode::TT_ENZYMATIC || _mode == STRMode::TP_ENZYMATIC;
    }

    const bool is_fusion() const {
        return _mode == STRMode::TT_FUSION || _mode == STRMode::TT_FUSION_ALLOW_SELF || _mode == STRMode::TP_FUSION;
    }

    const scalar rate() const {
        return _rate;
    }

    const scalar radius() const {
        return _radius;
    }

    const bool allow_self_connection() const {
        return _mode == STRMode::TT_FUSION_ALLOW_SELF;
    }

    const STRMode &mode() const {
        return _mode;
    }

private:

    friend class STRParser;

    SpatialTopologyReaction() = default;

    std::string _name;
    util::particle_type_pair _types;
    util::particle_type_pair _types_to;
    topology_type_pair _top_types;
    topology_type_pair _top_types_to;
    scalar _rate{0};
    scalar _radius{0};
    STRMode _mode{STRMode::TP_ENZYMATIC};
};

class STRParser {
public:

    explicit STRParser(const TopologyRegistry &registry) : _topology_registry(registry) {};

    /**
     *  Pass descriptor of form for Topology<->Topology
     * - fusion type:
     *     name: T1 (p1) + T2 (p2) -> T3 (p3--p4) [self=true|false]
     * - enzymatic type:
     *     name: T1 (p1) + T2 (p2) -> T3 (p3) + T4 (p4)
     *
     * and for Topology<->Particle
     * - fusion type:
     *     name: T1 (p1) + (p2) -> T2 (p3--p4)
     * - enzymatic type:
     *     name: T1 (p1) + (p2) -> T2 (p3) + (p4)
     *
     * @param descriptor the descriptor
     * @param rate the rate
     * @param radius the radius
     * @return the parsed reaction object
     */
    SpatialTopologyReaction parse(const std::string &descriptor, scalar rate, scalar radius) const;

private:
    std::reference_wrapper<const TopologyRegistry> _topology_registry;
};

NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
