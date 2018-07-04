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
 * @file GraphTopology.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/macros.h>

#include "Topology.h"
#include "graph/Graph.h"
#include "reactions/reactions.h"
#include "TopologyRegistry.h"
#include "Utils.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
class StateModel;
NAMESPACE_BEGIN(top)

class GraphTopology : public Topology {
public:

    using topology_graph = graph::Graph;
    using topology_reaction_rate = scalar;
    using topology_reaction_rates = std::vector<topology_reaction_rate>;
    using types_vec = std::vector<ParticleTypeId>;
    using vertex = graph::Vertex;

    /**
     * Creates a new graph topology. An internal graph object will be created with vertices corresponding to the
     * particles handed in.
     * @param type the type
     * @param particles the particles
     * @param types particle's types
     * @param context the kernel's context
     * @param stateModel the kernel's state model
     */
    GraphTopology(TopologyTypeId type, const particle_indices &particles, const types_vec &types,
                  const model::Context &context, const StateModel* stateModel);

    /**
     * Will create a graph topology out of an already existing graph and a list of particles, where the i-th vertex
     * of the graph will map to the i-th particles in the particles list.
     * @param type the type
     * @param particles the particles list
     * @param graph the already existing graph
     * @param context the kernel's context
     * @param stateModel the kernels state model
     */
    GraphTopology(TopologyTypeId type, particle_indices &&particles, topology_graph &&graph,
                  const model::Context &context, const StateModel* stateModel);

    ~GraphTopology() override = default;

    GraphTopology(GraphTopology &&) = default;

    GraphTopology &operator=(GraphTopology &&) = default;

    GraphTopology(const GraphTopology &) = delete;

    GraphTopology &operator=(const GraphTopology &) = delete;

    topology_graph &graph() {
        return graph_;
    }

    const topology_graph &graph() const {
        return graph_;
    }

    void configure();

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
            throw std::invalid_argument(fmt::format("The graph is not connected! (GEXF representation: {})",
                                                    util::to_gexf(graph())));
        }
    }

    std::vector<GraphTopology> connectedComponents();

    const bool isDeactivated() const {
        return deactivated;
    }

    void deactivate() {
        deactivated = true;
    }

    const bool isNormalParticle(const Kernel &k) const;

    const topology_reaction_rate cumulativeRate() const {
        return _cumulativeRate;
    }

    void appendParticle(particle_index newParticle, ParticleTypeId newParticleType,
                        particle_index counterPart, ParticleTypeId counterPartType);

    void appendTopology(GraphTopology &other, particle_index otherParticle, ParticleTypeId otherNewParticleType,
                        particle_index thisParticle, ParticleTypeId thisNewParticleType, TopologyTypeId newType);

    const TopologyTypeId &type() const {
        return _topology_type;
    }

    TopologyTypeId &type() {
        return _topology_type;
    }

    topology_graph::vertex_ref vertexForParticle(particle_index particle) {
        auto it = std::find(particles.begin(), particles.end(), particle);
        if(it != particles.end()) {
            return std::next(graph_.vertices().begin(), std::distance(particles.begin(), it));
        }
        return graph_.vertices().end();
    }

    const topology_reaction_rates &rates() const {
        return _reaction_rates;
    }

    const model::Context &context() const {
        return _context.get();
    }

    std::vector<Particle> fetchParticles() const;

    Particle particleForVertex(const vertex &vertex) const;

    Particle particleForVertex(topology_graph::vertex_ref vertexRef) const {
        return particleForVertex(*vertexRef);
    }

protected:
    topology_graph graph_;
    std::reference_wrapper<const model::Context> _context;
    const model::StateModel *_stateModel;
    topology_reaction_rates _reaction_rates;
    topology_reaction_rate _cumulativeRate;
    TopologyTypeId _topology_type;
    bool deactivated{false};
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
