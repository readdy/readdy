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

GraphTopology::GraphTopology(const Topology::particles_t &particles, const std::vector<particle_type_type> &types,
                             const graph::PotentialConfiguration *const config)
        : Topology(particles), config(config), graph_(std::make_unique<graph::Graph>()) {
    assert(types.size() == particles.size());
    std::size_t i = 0;
    for (auto itTypes = types.begin(); itTypes != types.end(); ++itTypes, ++i) {
        graph().addVertex(i, *itTypes);
    }
}

graph::Graph &GraphTopology::graph() {
    return *graph_;
}

const graph::Graph &GraphTopology::graph() const {
    return *graph_;
}

void GraphTopology::configure() {
    bondedPotentials.clear();
    anglePotentials.clear();
    torsionPotentials.clear();

    std::unordered_map<graph::BondType, std::vector<BondConfiguration>, readdy::util::hash::EnumClassHash> bonds;
    std::unordered_map<graph::AngleType, std::vector<AngleConfiguration>, readdy::util::hash::EnumClassHash> angles;
    std::unordered_map<graph::TorsionType, std::vector<DihedralConfiguration>, readdy::util::hash::EnumClassHash> dihedrals;

    graph_->findNTuples([&](const graph::Graph::vertex_ptr_tuple &tuple) {
        auto it = config->pairPotentials.find(
                std::tie(std::get<0>(tuple)->particleType(), std::get<1>(tuple)->particleType()));
        if (it != config->pairPotentials.end()) {
            for (const auto &cfg : it->second) {
                bonds[cfg.type].emplace_back(std::get<0>(tuple)->particleIndex, std::get<1>(tuple)->particleIndex,
                                             cfg.forceConstant, cfg.length);
            }
        }
    }, [&](const graph::Graph::vertex_ptr_triple &triple) {
        const auto &v1 = std::get<0>(triple);
        const auto &v2 = std::get<1>(triple);
        const auto &v3 = std::get<2>(triple);
        auto it = config->anglePotentials.find(std::tie(v1->particleType(), v2->particleType(), v3->particleType()));
        if (it != config->anglePotentials.end()) {
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
        auto it = config->torsionPotentials.find(
                std::tie(v1->particleType(), v2->particleType(), v3->particleType(), v4->particleType()));
        if (it != config->torsionPotentials.end()) {
            for (const auto &cfg : it->second) {
                dihedrals[cfg.type].emplace_back(v1->particleIndex, v2->particleIndex, v3->particleIndex,
                                                 v4->particleIndex, cfg.forceConstant, cfg.multiplicity,
                                                 cfg.phi_0);
            }
        }
    });
    for(const auto& bond : bonds) {
        switch (bond.first) {
            case graph::BondType::HARMONIC: {
                addBondedPotential(std::make_unique<HarmonicBondPotential>(this, bond.second));
                break;
            };
        }
    }
    for(const auto& angle : angles) {
        switch(angle.first) {
            case graph::AngleType::HARMONIC: {
                addAnglePotential(std::make_unique<HarmonicAnglePotential>(this, angle.second));
                break;
            };
        }
    }
    for(const auto& dih : dihedrals) {
        switch(dih.first) {
            case graph::TorsionType::COS_DIHEDRAL: {
                addTorsionPotential(std::make_unique<CosineDihedralPotential>(this, dih.second));
                break;
            };
        }
    }
}

void GraphTopology::validate() {
    if (!graph().isConnected()) {
        log::warn("The graph is not connected!");
    }
}

}
}
}

