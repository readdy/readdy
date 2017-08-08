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
    using graph_t = graph::Graph;

    using edge = graph_t::edge;
    using vertex = graph_t::vertex_ref;

    explicit TopologyReactionAction(GraphTopology *topology);

    virtual ~TopologyReactionAction() = default;

    TopologyReactionAction(const TopologyReactionAction&) = default;

    TopologyReactionAction& operator=(const TopologyReactionAction&) = delete;

    TopologyReactionAction(TopologyReactionAction&&) = default;

    TopologyReactionAction& operator=(TopologyReactionAction&&) = delete;

    virtual void execute() = 0;

    virtual void undo() = 0;

protected:
    GraphTopology *const topology;
};

class ChangeParticleType : public TopologyReactionAction {
public:

    ChangeParticleType(GraphTopology *topology, const vertex &v, const particle_type_type &type_to);

protected:
    vertex _vertex;
    particle_type_type type_to, previous_type;
};

class AddEdge : public TopologyReactionAction {
public:
    AddEdge(GraphTopology *topology, const edge &edge);

    void execute() override;

    void undo() override;

private:
    edge label_edge_;
};

class RemoveEdge : public TopologyReactionAction {
public:
    RemoveEdge(GraphTopology *topology, const edge &edge);

    void execute() override;

    void undo() override;

private:
    edge label_edge_;
};

NAMESPACE_END(actions)
NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
