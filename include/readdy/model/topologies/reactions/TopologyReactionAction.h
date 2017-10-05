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
 * @file TopologyReactionOperation.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <memory>
#include <readdy/common/macros.h>
#include <readdy/model/topologies/graph/Graph.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
class GraphTopology;
NAMESPACE_BEGIN(reactions)
NAMESPACE_BEGIN(actions)

class TopologyReactionAction {
public:
    /**
     * the GraphTopology's graph
     */
    using topology_graph = graph::Graph;

    /**
     * an edge in the graph
     */
    using edge = topology_graph::edge;
    /**
     * a vertex reference of the graph
     */
    using vertex = topology_graph::vertex_ref;

    /**
     * creates a new topology reaction action w.r.t. the given topology
     * @param topology the topology
     */
    explicit TopologyReactionAction(GraphTopology *topology);

    /**
     * default delete
     */
    virtual ~TopologyReactionAction() = default;

    /**
     * default copy
     */
    TopologyReactionAction(const TopologyReactionAction&) = default;

    /**
     * no copy assign
     */
    TopologyReactionAction& operator=(const TopologyReactionAction&) = delete;

    /**
     * default move
     */
    TopologyReactionAction(TopologyReactionAction&&) = default;

    /**
     * no move assign
     */
    TopologyReactionAction& operator=(TopologyReactionAction&&) = delete;

    /**
     * base method for executing the action
     */
    virtual void execute() = 0;

    /**
     * base method for undoing the action
     */
    virtual void undo() = 0;

protected:
    /**
     * a pointer to the topology on which this action should be executed
     */
    GraphTopology *const topology;
};

class ChangeParticleType : public TopologyReactionAction {
public:

    /**
     * Creates an action that changes the particle type of one of the particles in the topology. Since the data
     * structures may vary with the selected kernel, the implementation of do/undo is kernel specific.
     *
     * @param topology the respective topology
     * @param v the vertex pointing to the particle whose type should be changed
     * @param type_to the target type
     */
    ChangeParticleType(GraphTopology *topology, const vertex &v, const particle_type_type &type_to);

protected:
    /**
     * a reference to the vertex
     */
    vertex _vertex;
    /**
     * the target type
     */
    particle_type_type type_to;
    /**
     * the previous particle type, stored for undo
     */
    particle_type_type previous_type;
};

class AddEdge : public TopologyReactionAction {
public:
    /**
     * Creates an action that introduces an edge in the graph.
     * @param topology the topology
     * @param edge the edge to introduce
     */
    AddEdge(GraphTopology *topology, const edge &edge);

    /**
     * do!
     */
    void execute() override;

    /**
     * undo!
     */
    void undo() override;

private:
    /**
     * the edge to introduce
     */
    edge label_edge_;
};

class RemoveEdge : public TopologyReactionAction {
public:
    /**
     * Creates an action that removes an edge in the topology.
     * @param topology the topology
     * @param edge the edge to remove
     */
    RemoveEdge(GraphTopology *topology, const edge &edge);

    /**
     * execute me
     */
    void execute() override;

    /**
     * oops, undo
     */
    void undo() override;

private:
    /**
     * the edge to remove
     */
    edge label_edge_;
};

class ChangeTopologyType : public TopologyReactionAction {
public:
    /**
     * Creates an action that changes the topology type of this topology
     * @param topology this topology
     * @param newType the target type
     */
    ChangeTopologyType(GraphTopology *topology, topology_type_type newType);

    /**
     * execute me
     */
    void execute() override;

    /**
     * oops, undo
     */
    void undo() override;

private:
    /**
     * the target type
     */
    topology_type_type _newType;
    /**
     * the previous type, stored for undo
     */
    topology_type_type _prevType {0};
};

NAMESPACE_END(actions)
NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
