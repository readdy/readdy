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
 * @file GraphTopology.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.03.17
 * @copyright BSD-3
 */

#pragma once

#include <graphs/graphs.h>

#include "common.h"
#include "Topology.h"
#include "reactions/reactions.h"
#include "TopologyRegistry.h"

namespace readdy::model {
class StateModel;
namespace top {

class GraphTopology : public Topology {
public:
    using Vertex = graphs::Vertex<VertexData>;
    using Graph = graphs::Graph<Vertex>;
    using ReactionRate = scalar;
    using ReactionRates = std::vector<ReactionRate>;

    /**
     * Creates a new graph topology. An internal graph object will be created with vertices corresponding to the
     * particles handed in.
     * @param type the type
     * @param context the kernel's context
     * @param stateModel the kernel's state model
     */
    GraphTopology(TopologyTypeId type, Graph graph, const model::Context &context, const StateModel *stateModel);

    ~GraphTopology() override = default;

    GraphTopology(GraphTopology &&) = default;

    GraphTopology &operator=(GraphTopology &&) = default;

    GraphTopology(const GraphTopology &) = delete;

    GraphTopology &operator=(const GraphTopology &) = delete;

    void setGraph(Graph graph) {
        _graph = std::move(graph);
    }

    [[nodiscard]] const Graph &graph() const {
        return _graph;
    }

    [[nodiscard]] auto containsEdge(Graph::Edge edge) const {
        return _graph.containsEdge(std::forward<Graph::Edge>(edge));
    }

    [[nodiscard]] auto containsEdge(Graph::VertexIndex ix1, Graph::VertexIndex ix2) const {
        return _graph.containsEdge(ix1, ix2);
    }

    void addEdge(Graph::Edge edge) {
        _graph.addEdge(edge);
    }

    void addEdge(Graph::VertexIndex ix1, Graph::VertexIndex ix2) {
        _graph.addEdge(ix1, ix2);
    }

    void removeEdge(Graph::Edge edge) {
        _graph.removeEdge(edge);
    }

    void removeEdge(Graph::VertexIndex ix1, Graph::VertexIndex ix2) {
        _graph.removeEdge(ix1, ix2);
    }

    void configure();

    [[nodiscard]] ParticleTypeId typeOf(Graph::VertexIndex vertex) const;

    [[nodiscard]] ParticleTypeId typeOf(const Vertex &v) const;

    void updateReactionRates(const TopologyRegistry::StructuralReactionCollection &reactions) {
        _cumulativeRate = 0;
        _reaction_rates.resize(reactions.size());
        auto that = this;
        std::transform(reactions.begin(), reactions.end(), _reaction_rates.begin(), [&](const auto &reaction) {
            const auto rate = reaction.rate(*that);
            _cumulativeRate += rate;
            return rate;
        });
    }

    void validate() {
        if (!graph().isConnected()) {
            throw std::invalid_argument(fmt::format("The graph is not connected! (GEXF representation: {})", _graph.gexf()));
        }
        for(const auto [i1, i2] : graph().edges()) {
            if(graph().vertices().at(i1).deactivated()) {
                throw std::logic_error(fmt::format("Edge ({} -- {}) points to deactivated vertex {}!", i1, i2, i1));
            }
            if(graph().vertices().at(i2).deactivated()) {
                throw std::logic_error(fmt::format("Edge ({} -- {}) points to deactivated vertex {}!", i1, i2, i2));
            }
        }
        for(const auto& p : fetchParticles()) {
            if(context().particleTypes().infoOf(p.type()).flavor != particleflavor::TOPOLOGY) {
                throw std::runtime_error(fmt::format("Topology contains particle {} which is not a topology particle!", p));
            }
        }
    }

    std::vector<GraphTopology> connectedComponents();

    [[nodiscard]] bool isDeactivated() const {
        return deactivated;
    }

    void deactivate() {
        deactivated = true;
    }

    [[nodiscard]] bool isNormalParticle(const Kernel &k) const;

    [[nodiscard]] ReactionRate cumulativeRate() const {
        return _cumulativeRate;
    }

    [[nodiscard]] typename Graph::VertexList::iterator vertexIteratorForParticle(VertexData::ParticleIndex index);

    [[nodiscard]] typename Graph::VertexList::const_iterator vertexIteratorForParticle(VertexData::ParticleIndex index) const;

    [[nodiscard]] Graph::VertexIndex vertexIndexForParticle(VertexData::ParticleIndex index) const {
        auto it = vertexIteratorForParticle(index);
        if(it == _graph.vertices().end()) {
            throw std::invalid_argument(fmt::format("Particle {} not contained in graph", index));
        }
        return std::distance(_graph.vertices().begin(), it);
    }

    typename Graph::VertexIndex appendParticle(VertexData::ParticleIndex newParticle,
                                               VertexData::ParticleIndex counterPart);

    typename Graph::VertexIndex appendTopology(const GraphTopology &other, VertexData::ParticleIndex otherParticle,
                                               VertexData::ParticleIndex thisParticle, TopologyTypeId newType);

    [[nodiscard]] std::vector<VertexData::ParticleIndex> particleIndices() const;

    [[nodiscard]] TopologyTypeId type() const {
        return _topology_type;
    }

    TopologyTypeId &type() {
        return _topology_type;
    }

    [[nodiscard]] const ReactionRates &rates() const {
        return _reaction_rates;
    }

    [[nodiscard]] const model::Context &context() const {
        return _context.get();
    }

    [[nodiscard]] std::vector<Particle> fetchParticles() const;

    [[nodiscard]] Particle particleForVertex(const Vertex &vertex) const;

    [[nodiscard]] Particle particleForVertex(Graph::VertexIndex vertexRef) const {
        return particleForVertex(_graph.vertices().at(vertexRef));
    }

    [[nodiscard]] auto nParticles() const {
        return _graph.nVertices();
    }

protected:
    Graph _graph;
    std::reference_wrapper<const model::Context> _context;
    const model::StateModel *_stateModel;
    ReactionRates _reaction_rates;
    ReactionRate _cumulativeRate{};
    TopologyTypeId _topology_type;
    bool deactivated{false};
};

}
}
