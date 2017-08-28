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

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

class GraphTopology : public Topology {
public:

    using topology_graph = graph::Graph;
    using topology_reaction_rate = scalar;
    using topology_reaction_rates = std::vector<topology_reaction_rate>;
    using types_vec = std::vector<particle_type_type>;

    /**
     * Creates a new graph topology. An internal graph object will be created with vertices corresponding to the
     * particles handed in.
     * @param type the type
     * @param particles the particles
     * @param types particle's types
     * @param config the configuration table
     */
    GraphTopology(topology_type_type type, const particle_indices &particles, const types_vec &types,
                  const api::PotentialConfiguration &config);

    /**
     * Will create a graph topology out of an already existing graph and a list of particles, where the i-th vertex
     * of the graph will map to the i-th particles in the particles list.
     * @param type the type
     * @param particles the particles list
     * @param graph the already existing graph
     * @param config the configuration table
     */
    GraphTopology(topology_type_type type, particle_indices &&particles, topology_graph &&graph,
                  const api::PotentialConfiguration &config);

    virtual ~GraphTopology() = default;

    GraphTopology(GraphTopology &&) = default;

    GraphTopology &operator=(GraphTopology &&) = default;

    GraphTopology(const GraphTopology &) = delete;

    GraphTopology &operator=(const GraphTopology &) = delete;

    topology_graph &graph();

    const topology_graph &graph() const;

    void configure();

    void updateReactionRates(const TopologyRegistry::structural_reactions &reactions);

    void validate();

    std::vector<GraphTopology> connectedComponents();

    const bool isDeactivated() const;

    void deactivate();

    const bool isNormalParticle(const Kernel &k) const;

    const topology_reaction_rate cumulativeRate() const;

    void appendParticle(particle_index newParticle, particle_type_type newParticleType,
                        particle_index counterPart, particle_type_type counterPartType);

    void appendTopology(GraphTopology &other, particle_index otherParticle, particle_type_type otherNewParticleType,
                        particle_index thisParticle, particle_type_type thisNewParticleType);

    topology_type_type type() const;

    topology_graph::vertex_ref vertexForParticle(particle_index particle);

    const topology_reaction_rates &rates() const;

protected:
    topology_graph graph_;
    std::reference_wrapper<const api::PotentialConfiguration> config;
    topology_reaction_rates _reaction_rates;
    topology_reaction_rate _cumulativeRate;
    topology_type_type _topology_type;
    bool deactivated{false};
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
