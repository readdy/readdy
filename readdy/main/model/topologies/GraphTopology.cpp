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
#include <sstream>


namespace readdy {
namespace model {
namespace top {

GraphTopology::GraphTopology(const Topology::particles_t &particles, const std::vector<particle_type_type> &types,
                             std::reference_wrapper<const api::PotentialConfiguration> config)
        : Topology(particles), config(config), graph_() {
    assert(types.size() == particles.size());
    std::size_t i = 0;
    for (auto itTypes = types.begin(); itTypes != types.end(); ++itTypes, ++i) {
        graph().addVertex(i, *itTypes);
    }
}

GraphTopology::GraphTopology(const Topology::particles_t &particles, graph::Graph &&graph,
                             std::reference_wrapper<const api::PotentialConfiguration> config)
        : Topology(particles), config(config), graph_(std::move(graph)) {
    assert(GraphTopology::graph().vertices().size() == particles.size());
    if (GraphTopology::graph().vertices().size() != particles.size()) {
        throw std::invalid_argument("the number of particles and the number of vertices should match when creating"
                                            "a graph in this way!");
    }
    std::size_t idx = 0;
    auto &vertices = GraphTopology::graph().vertices();
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        it->particleIndex = idx++;
    }
}

graph::Graph &GraphTopology::graph() {
    return graph_;
}

const graph::Graph &GraphTopology::graph() const {
    return graph_;
}

void GraphTopology::configure() {

    validate();

    bondedPotentials.clear();
    anglePotentials.clear();
    torsionPotentials.clear();

    std::unordered_map<api::BondType, std::vector<BondConfiguration>, readdy::util::hash::EnumClassHash> bonds;
    std::unordered_map<api::AngleType, std::vector<AngleConfiguration>, readdy::util::hash::EnumClassHash> angles;
    std::unordered_map<api::TorsionType, std::vector<DihedralConfiguration>, readdy::util::hash::EnumClassHash> dihedrals;

    graph_.findNTuples([&](const graph::Graph::vertex_ptr_tuple &tuple) {
        auto v1 = std::get<0>(tuple);
        auto v2 = std::get<1>(tuple);
        auto it = config.pairPotentials.find(
                std::tie(v1->particleType(), v2->particleType()));
        if (it != config.pairPotentials.end()) {
            for (const auto &cfg : it->second) {
                bonds[cfg.type].emplace_back(v1->particleIndex, v2->particleIndex,
                                             cfg.forceConstant, cfg.length);
            }
        } else {
            std::ostringstream ss;

            ss << "The edge " << v1->particleIndex;
            if (!v1->label.empty()) {
                ss << " (" << v1->label << ")";
            }
            ss << " -- " << v2->particleIndex;
            if (!v2->label.empty()) {
                ss << " (" << v2->label << ")";
            }
            ss << " has no bond configured! (See KernelContext.configureTopologyBondPotential())";

            throw std::invalid_argument(ss.str());
        }
    }, [&](const graph::Graph::vertex_ptr_triple &triple) {
        const auto &v1 = std::get<0>(triple);
        const auto &v2 = std::get<1>(triple);
        const auto &v3 = std::get<2>(triple);
        auto it = config.anglePotentials.find(std::tie(v1->particleType(), v2->particleType(), v3->particleType()));
        if (it != config.anglePotentials.end()) {
            for (const auto &cfg : it->second) {
                angles[cfg.type].emplace_back(v1->particleIndex, v2->particleIndex, v3->particleIndex,
                                              cfg.forceConstant, cfg.equilibriumAngle);
            }
        }
    }, [&](const graph::Graph::vertex_ptr_quadruple &quadruple) {
        const auto &v1 = std::get<0>(quadruple);
        const auto &v2 = std::get<1>(quadruple);
        const auto &v3 = std::get<2>(quadruple);
        const auto &v4 = std::get<3>(quadruple);
        auto it = config.torsionPotentials.find(
                std::tie(v1->particleType(), v2->particleType(), v3->particleType(), v4->particleType()));
        if (it != config.torsionPotentials.end()) {
            for (const auto &cfg : it->second) {
                dihedrals[cfg.type].emplace_back(v1->particleIndex, v2->particleIndex, v3->particleIndex,
                                                 v4->particleIndex, cfg.forceConstant, cfg.multiplicity,
                                                 cfg.phi_0);
            }
        }
    });
    for (const auto &bond : bonds) {
        switch (bond.first) {
            case api::BondType::HARMONIC: {
                addBondedPotential(std::make_unique<HarmonicBondPotential>(this, bond.second));
                break;
            };
        }
    }
    for (const auto &angle : angles) {
        switch (angle.first) {
            case api::AngleType::HARMONIC: {
                addAnglePotential(std::make_unique<HarmonicAnglePotential>(this, angle.second));
                break;
            };
        }
    }
    for (const auto &dih : dihedrals) {
        switch (dih.first) {
            case api::TorsionType::COS_DIHEDRAL: {
                addTorsionPotential(std::make_unique<CosineDihedralPotential>(this, dih.second));
                break;
            };
        }
    }
}

void GraphTopology::validate() {
    if (!graph().isConnected()) {
        throw std::invalid_argument("The graph is not connected!");
    }
}

}
}
}

