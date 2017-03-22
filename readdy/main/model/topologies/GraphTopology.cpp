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
    auto itParticles = particles.begin();
    auto itTypes = types.begin();
    for (; itParticles != particles.end(); ++itParticles, ++itTypes) {
        graph().addVertex(*itParticles, *itTypes);
    }
}

graph::Graph &GraphTopology::graph() {
    return *graph_;
}

const graph::Graph &GraphTopology::graph() const {
    return *graph_;
}

void GraphTopology::configure() {
    std::unordered_map<graph::BondType, std::vector<BondConfiguration>> bonds;
    std::unordered_map<graph::AngleType, std::vector<AngleConfiguration>> angles;
    std::unordered_map<graph::TorsionType, std::vector<DihedralConfiguration>> dihedrals;
    for (auto &v : graph().vertices()) {
        v.visited = false;
    }

    for (auto &v : graph().vertices()) {
        v.visited = true;
        auto v_type = v.particleType();
        auto v_idx = v.particleIndex;
        // todo this is wrong
        for (auto vv : v.neighbors()) {
            if (!vv->visited) {
                vv->visited = true;
                auto vv_type = vv->particleType();
                auto vv_idx = vv->particleIndex;
                // got tuple v, vv
                {
                    log::debug("got type tuple ({}, {}) for particles {}, {}", v_type, vv_type, v_idx, vv_idx);
                    auto it = config->pairPotentials.find(std::tie(v_type, vv_type));
                    if (it != config->pairPotentials.end()) {
                        for (const auto &cfg : it->second) {
                            bonds[cfg.type].emplace_back(v_idx, vv_idx, cfg.forceConstant, cfg.length);
                        }
                    }
                }
                for (auto vvv : vv->neighbors()) {
                    if (!vvv->visited) {
                        vvv->visited = true;
                        auto vvv_type = vvv->particleType();
                        auto vvv_idx = vvv->particleIndex;
                        // got triple v, vv, vvv
                        log::debug("got type triple ({}, {}, {}) for particles {}, {}, {}", v_type, vv_type, vvv_type,
                                   v_idx, vv_idx, vvv_idx);
                        {
                            auto it = config->anglePotentials.find(std::tie(v_type, vv_type, vvv_type));
                            if (it != config->anglePotentials.end()) {
                                for (const auto &cfg : it->second) {
                                    angles[cfg.type].emplace_back(v_idx, vv_idx, vvv_idx, cfg.forceConstant,
                                                                  cfg.equilibriumAngle);
                                }
                            }
                        }
                        for (auto vvvv : vvv->neighbors()) {
                            if (!vvvv->visited) {
                                auto vvvv_type = vvvv->particleType();
                                auto vvvv_idx = vvvv->particleIndex;
                                // got quadruple v, vv, vvv, vvvv
                                log::debug("got type quadruple ({}, {}, {}, {}) for particles {}, {}, {}, {}", v_type,
                                           vv_type, vvv_type, vvvv_type, v_idx, vv_idx, vvv_idx, vvvv_idx);
                                {
                                    auto it = config->torsionPotentials.find(
                                            std::tie(v_type, vv_type, vvv_type, vvvv_type));
                                    if (it != config->torsionPotentials.end()) {
                                        for (const auto &cfg : it->second) {
                                            dihedrals[cfg.type].emplace_back(v_idx, vv_idx, vvv_idx, vvvv_idx,
                                                                             cfg.forceConstant, cfg.multiplicity,
                                                                             cfg.phi_0);
                                        }
                                    }
                                }
                            }
                        }
                        vvv->visited = false;
                    }
                }
                vv->visited = false;
            }
        }
    }
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

void GraphTopology::permuteIndices(const std::vector<std::size_t> &permutation) {
    Topology::permuteIndices(permutation);
    for (auto &vertex : graph().vertices()) {
        vertex.particleIndex = permutation[vertex.particleIndex];
    }
}

}
}
}

