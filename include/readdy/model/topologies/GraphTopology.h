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
#include "connectivity/Graph.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

class GraphTopology : public Topology {
public:

    using vertex_ptr = graph::Graph::vertices_t::iterator;
    using vertex_ptr_tuple = std::tuple<vertex_ptr, vertex_ptr>;
    using vertex_ptr_triple = std::tuple<vertex_ptr, vertex_ptr, vertex_ptr>;
    using vertex_ptr_quadruple = std::tuple<vertex_ptr, vertex_ptr, vertex_ptr, vertex_ptr>;

    GraphTopology(const particles_t& particles, const std::vector<particle_type_type> &types,
                  const graph::PotentialConfiguration *const config);

    virtual ~GraphTopology() = default;

    GraphTopology(GraphTopology &&) = delete;

    GraphTopology &operator=(GraphTopology &&) = delete;

    GraphTopology(const GraphTopology &) = delete;

    GraphTopology &operator=(const GraphTopology &) = delete;

    graph::Graph &graph();

    const graph::Graph &graph() const;

    void findNTuples(const std::function<void(const vertex_ptr_tuple&)> &tuple_callback,
                     const std::function<void(const vertex_ptr_triple&)> &triple_callback,
                     const std::function<void(const vertex_ptr_quadruple&)> &quadruple_callback);

    void configure();

    void validate();

    virtual void permuteIndices(const std::vector<std::size_t> &permutation) override;

private:
    std::unique_ptr<graph::Graph> graph_;
    const graph::PotentialConfiguration *const config;
};


NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
