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
#include <readdy/model/topologies/graph/Vertex.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
class GraphTopology;
NAMESPACE_BEGIN(reactions)
NAMESPACE_BEGIN(op)

class Operation;

using OperationRef = std::shared_ptr<Operation>;

class Operation {
public:
    Operation(GraphTopology *const topology);

    virtual ~Operation() = default;

    virtual void execute() = 0;

    virtual void undo() = 0;

protected:
    GraphTopology *const topology;
};

class ChangeParticleType : public Operation {
public:

    ChangeParticleType(GraphTopology *const topology, const graph::vertex_ref &v, const particle_type_type &type_to);

    void execute() override;

    void undo() override;

private:
    graph::vertex_ref vertex;
    particle_type_type type_to, previous_type;
};

class AddEdge : public Operation {
public:
    AddEdge(GraphTopology *const topology, const graph::edge &edge);

    void execute() override;

    void undo() override;

private:
    graph::edge edge;
};

class RemoveEdge : public Operation {
public:
    RemoveEdge(GraphTopology *const topology, const graph::edge &edge);

    void execute() override;

    void undo() override;

private:
    graph::edge edge;
};

NAMESPACE_END(op)
NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)