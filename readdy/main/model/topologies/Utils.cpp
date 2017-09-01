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
 * @file Utils.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/Utils.h>
#include <sstream>

namespace readdy {
namespace model {
namespace top {
namespace util {

std::string to_gexf(graph::Graph& graph) {
    std::ostringstream ss;
    ss << R"(<?xml version="1.0" encoding="UTF-8"?>)";
    ss << R"(<gexf xmlns="http://www.gexf.net/1.2draft" version="1.2">)";
    ss << R"(<graph mode="static" defaultedgetype="undirected">)";
    {
        ss << "<nodes>";
        std::size_t id = 0;
        for (auto &v : graph.vertices()) {
            ss << "<node id=\"" << v.particleIndex << "\"/>";
            v.visited = false;
            ++id;
        }
        ss << "</nodes>";
    }
    {
        ss << "<edges>";
        std::size_t id = 0;
        for(auto& v : graph.vertices()) {
            for(const auto& neighbor : v.neighbors()) {
                if(!neighbor->visited) {
                    ss << "<edge id=\"" << id << "\" "
                            "source=\"" << v.particleIndex << "\" "
                            "target=\"" << neighbor->particleIndex << "\" />";
                    ++id;
                }
            }
            v.visited = true;
        }
        ss << "</edges>";
    }
    ss << "</graph>";
    ss << "</gexf>";
    return ss.str();
}

}
}
}
}