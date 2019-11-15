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
 * << detailed description >>
 *
 * @file Operations.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 05.04.17
 * @copyright BSD-3
 */

#include <readdy/model/topologies/reactions/TopologyReactionAction.h>
#include <readdy/model/topologies/GraphTopology.h>

#include <utility>

namespace readdy::model::top::reactions::actions {

TopologyReactionAction::TopologyReactionAction(GraphTopology *const topology) : topology(topology){ }

ChangeParticleType::ChangeParticleType(GraphTopology *const topology, Graph::VertexIndex v,
                                       ParticleTypeId type_to)
        : TopologyReactionAction(topology), _vertex(v), type_to(type_to), previous_type(type_to){}

AddEdge::AddEdge(GraphTopology *const topology, Graph::Edge edge)
        : TopologyReactionAction(topology), label_edge_(std::move(edge)) {}

void AddEdge::execute() {
    topology->addEdge(label_edge_);
}

RemoveEdge::RemoveEdge(GraphTopology *const topology, Graph::Edge edge)
        : TopologyReactionAction(topology), label_edge_(std::move(edge)) {}

void RemoveEdge::execute() {
    topology->removeEdge(label_edge_);
}


ChangeTopologyType::ChangeTopologyType(GraphTopology *topology, TopologyTypeId newType)
        : TopologyReactionAction(topology), _newType(newType){}

void ChangeTopologyType::execute() {
    topology->type() = _newType;
}

ChangeParticlePosition::ChangeParticlePosition(GraphTopology *topology, Graph::VertexIndex v, Vec3 posTo)
        : TopologyReactionAction(topology), _vertex(v), _posTo(posTo) {}

AppendParticle::AppendParticle(GraphTopology *topology, std::vector<Graph::VertexIndex> neighbors,
                               ParticleTypeId type, Vec3 pos)
        : TopologyReactionAction(topology), neighbors(std::move(neighbors)), type(type), pos(pos) {}
}
