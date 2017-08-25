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
 * @file TopologyTypeRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 24.08.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/common.h>
#include <readdy/model/topologies/reactions/TopologyReaction.h>


NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

struct TopologyTypeInfo {
    using internal_topology_reaction = reactions::TopologyReaction;
    using internal_topology_reactions = std::vector<internal_topology_reaction>;

    std::string name {""};
    topology_type_type type {0};
    internal_topology_reactions internal_reactions {};
};

class TopologyTypeRegistry {
public:

    using internal_topology_reaction = TopologyTypeInfo::internal_topology_reaction;
    using internal_topology_reactions = TopologyTypeInfo::internal_topology_reactions;

    using topology_type_registry = std::unordered_map<topology_type_type, TopologyTypeInfo>;

    TopologyTypeRegistry() = default;
    TopologyTypeRegistry(const TopologyTypeRegistry&) = delete;
    TopologyTypeRegistry& operator=(const TopologyTypeRegistry&) = delete;
    TopologyTypeRegistry(TopologyTypeRegistry&&) = delete;
    TopologyTypeRegistry& operator=(TopologyTypeRegistry&&) = delete;
    ~TopologyTypeRegistry() = default;

    topology_type_type add_type(const std::string& name, const internal_topology_reactions& reactions = {});

    void add_reaction(topology_type_type type, const reactions::TopologyReaction &reaction);

    void add_reaction(topology_type_type type, reactions::TopologyReaction &&reaction);

    void add_reaction(const std::string &type, const reactions::TopologyReaction &reaction);

    void add_reaction(const std::string &type, reactions::TopologyReaction &&reaction);

    const TopologyTypeInfo &info_of(topology_type_type type) const;

    const std::string &name_of(topology_type_type type) const;

    topology_type_type id_of(const std::string& name) const;

    bool empty();

    bool contains_reactions();

    const TopologyTypeInfo::internal_topology_reactions &reactions_of(topology_type_type type) const;

    const internal_topology_reactions &reactions_of(const std::string &type) const;

    void configure();

    void debug_output() const;

private:
    static unsigned short counter;

    TopologyTypeInfo defaultInfo;
    bool _contains_reactions {false};
    topology_type_registry _registry {};

};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)