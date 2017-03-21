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
    GraphTopology(particles_t &&, const graph::PotentialConfiguration *const config);

    GraphTopology(GraphTopology&&) = delete;
    GraphTopology& operator=(GraphTopology&&) = delete;
    GraphTopology(const GraphTopology&) = delete;
    GraphTopology& operator=(const GraphTopology&) = delete;

    graph::Graph* const graph();

    const graph::Graph* const graph() const;

    void configureByGraph();

private:
    std::unique_ptr<graph::Graph> graph_;
    const graph::PotentialConfiguration* const config;
};


NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
