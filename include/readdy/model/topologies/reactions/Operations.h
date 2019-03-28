/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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
 * @copyright BSD-3
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
    ChangeParticleType(const vertex_ref &vertex, ParticleTypeId type_to) : _vertex(vertex), _type_to(type_to) {}

    /**
     * Create the corresponding action.
     * @param topology the topology
     * @param factory the action factory
     * @return a pointer to the change particle type action
     */
    action_ptr create_action(topology_ref topology, factory_ref factory) const override {
        return factory->createChangeParticleType(topology, _vertex, _type_to);
    }

private:
    vertex_ref _vertex;
    ParticleTypeId _type_to;
};

class AppendParticle : public Operation {
public:

    AppendParticle(std::vector<vertex_ref> neighbors, ParticleTypeId type, const Vec3 &pos)
        : neighbors(std::move(neighbors)), type(type), pos(pos) {};

    action_ptr create_action(topology_ref topology, factory_ref factory) const override {
        return factory->createAppendParticle(topology, neighbors, type, pos);
    }

private:
    std::vector<vertex_ref> neighbors;
    ParticleTypeId type;
    Vec3 pos;
};

class ChangeParticlePosition : public Operation {
public:
    ChangeParticlePosition(const vertex_ref &vertex, Vec3 position) : _vertex(vertex), _pos(position) {};

    action_ptr create_action(topology_ref topology, factory_ref factory) const override {
        return factory->createChangeParticlePosition(topology, _vertex, _pos);
    }

private:
    vertex_ref _vertex;
    Vec3 _pos;
};

class ChangeTopologyType : public Operation {
public:
    /**
     * Creates an action that changes the topology type of the belonging topology.
     * @param type_to the target type
     */
    explicit ChangeTopologyType(std::string type_to) : _type_to(std::move(type_to)) {};

    /**
     * Create the corresponding action.
     * @param topology the topology
     * @param factory the action factory
     * @return a pointer to the action
     */
    action_ptr create_action(topology_ref topology, factory_ref factory) const override {
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
    explicit AddEdge(edge edge) : _edge(std::move(edge)) {};

    /**
     * Create the corresponding action
     * @param topology the topology
     * @param factory the action factory
     * @return a pointer to the respective action
     */
    action_ptr create_action(topology_ref topology, factory_ref factory) const override {
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
    explicit RemoveEdge(edge edge) : _edge(std::move(edge)) {};

    /**
     * Create the corresponding action
     * @param topology the topology
     * @param factory the action factory
     * @return a pointer to the respective action
     */
    action_ptr create_action(topology_ref topology, factory_ref factory) const override {
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
