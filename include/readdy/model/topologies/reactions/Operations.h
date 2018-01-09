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
 * This header files contains a base class for all operations on graph topologies, as well as definitions for:
 *   * ChangeParticleType
 *   * ChangeTopologyType
 *   * AddEdge
 *   * RemoveEdge
 *
 * @file Operations.h
 * @brief Definitions for various operations that can be performed on graph topologies.
 * @author clonker
 * @date 13.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <memory>
#include <readdy/common/macros.h>
#include "TopologyReactionAction.h"
#include "TopologyReactionActionFactory.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
class GraphTopology;
NAMESPACE_BEGIN(reactions)
NAMESPACE_BEGIN(op)

class Operation {
public:
    /**
     * Reference to an operation
     */
    using Ref = std::shared_ptr<Operation>;
    /**
     * Reference to the respective topology reaction action factory
     */
    using factory_ref = const actions::TopologyReactionActionFactory *const;
    /**
     * Reference to the respective graph topology
     */
    using topology_ref = GraphTopology *const;
    /**
     * Type of the graph of topologies
     */
    using topology_graph = actions::TopologyReactionAction::topology_graph;
    /**
     * pointer type to a topology reaction action
     */
    using action_ptr = std::unique_ptr<actions::TopologyReactionAction>;
    /**
     * reference to a vertex
     */
    using vertex_ref = topology_graph::vertex_ref;
    /**
     * an edge
     */
    using edge = topology_graph::edge;

    /**
     * Interface of the create_action method which will create the corresponding action on the selected kernel.
     * @param topology the topology this action should act upon
     * @param factory the factory
     * @return a unique pointer to the action
     */
    virtual action_ptr create_action(topology_ref topology, factory_ref factory) const = 0;
};

class ChangeParticleType : public Operation {
public:
    /**
     * Creates an action that changes the particle type of the particle pointed to by vertex.
     * @param vertex the vertex
     * @param type_to the target type
     */
    ChangeParticleType(const vertex_ref &vertex, particle_type_type type_to) : _vertex(vertex), _type_to(type_to) {}

    /**
     * Create the corresponding action.
     * @param topology the topology
     * @param factory the action factory
     * @return a pointer to the change particle type action
     */
    virtual action_ptr create_action(topology_ref topology, factory_ref factory) const override {
        return factory->createChangeParticleType(topology, _vertex, _type_to);
    }

private:
    vertex_ref _vertex;
    particle_type_type _type_to;
};

class ChangeTopologyType : public Operation {
public:
    /**
     * Creates an action that changes the topology type of the belonging topology.
     * @param type_to the target type
     */
    explicit ChangeTopologyType(const std::string &type_to) : _type_to(type_to) {};

    /**
     * Create the corresponding action.
     * @param topology the topology
     * @param factory the action factory
     * @return a pointer to the action
     */
    virtual action_ptr create_action(topology_ref topology, factory_ref factory) const override {
        return factory->createChangeTopologyType(topology, _type_to);
    }

private:
    std::string _type_to;
};

class AddEdge : public Operation {
public:
    /**
     * Adds the specified edge on the graph.
     * @param edge the edge
     */
    explicit AddEdge(const edge &edge) : _edge(edge) {};

    /**
     * Create the corresponding action
     * @param topology the topology
     * @param factory the action factory
     * @return a pointer to the respective action
     */
    virtual action_ptr create_action(topology_ref topology, factory_ref factory) const override {
        return factory->createAddEdge(topology, _edge);
    }

private:
    edge _edge;
};

class RemoveEdge : public Operation {
public:
    /**
     * Operation for removing the specified edge on the graph.
     * @param edge the edge
     */
    explicit RemoveEdge(const edge &edge) : _edge(edge) {};

    /**
     * Create the corresponding action
     * @param topology the topology
     * @param factory the action factory
     * @return a pointer to the respective action
     */
    virtual action_ptr create_action(topology_ref topology, factory_ref factory) const override {
        return factory->createRemoveEdge(topology, _edge);
    }

private:
    edge _edge;
};

NAMESPACE_END(op)
NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
