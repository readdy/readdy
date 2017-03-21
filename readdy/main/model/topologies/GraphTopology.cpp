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
 * @file GraphTopology.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/GraphTopology.h>


namespace readdy {
namespace model {
namespace top {

GraphTopology::GraphTopology(Topology::particles_t &&particles, const graph::PotentialConfiguration *const config)
        : Topology(std::move(particles)), config(config), graph_(std::make_unique<graph::Graph>()) {
    std::for_each(this->particles.begin(), this->particles.end(),
                  [this](std::size_t id) { graph()->addVertex({id}); });
}

graph::Graph *const GraphTopology::graph() {
    return graph_.get();
}

const graph::Graph *const GraphTopology::graph() const {
    return graph_.get();
}

void GraphTopology::configureByGraph() {
    if (!graph()) {
        log::critical("This should not be called if the topology was requested without graph!");
    } else {

    }
}

}
}
}

