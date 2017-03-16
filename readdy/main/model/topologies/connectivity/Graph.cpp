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
 * @file Graph.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 16.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/connectivity/Graph.h>

namespace readdy {
namespace model {
namespace top {
namespace graph {

void Graph::addEdge(std::size_t v1, std::size_t v2) {
    if (v1 > vertices.size() || v2 > vertices.size() || v1 == v2) {
        throw std::invalid_argument(
                "When adding an edge the start and end vertices must not be equal and the "
                        "indices be contained in the vertices vector. v1=" + std::to_string(v1) + ", v2=" +
                std::to_string(v2) + ", #vertices=" + std::to_string(vertices.size()));
    }
    auto& vertex1 = vertices.at(v1);
    auto& vertex2 = vertices.at(v2);
    vertex1.neighbors.push_back(&vertex2);
    vertex2.neighbors.push_back(&vertex1);
}


}
}
}
}