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

#include <sstream>

#include <readdy/model/Kernel.h>
#include <readdy/model/topologies/GraphTopology.h>


namespace readdy {
namespace model {
namespace top {

GraphTopology::GraphTopology(const Topology::particles_t &particles, const std::vector<particle_type_type> &types,
                             const api::PotentialConfiguration& config)
        : Topology(particles), config(config), graph_() {
    assert(types.size() == particles.size());
    std::size_t i = 0;
    for (auto itTypes = types.begin(); itTypes != types.end(); ++itTypes, ++i) {
        graph().addVertex(i, *itTypes);
    }
}

GraphTopology::GraphTopology(Topology::particles_t &&particles, graph::Graph &&graph,
                             const api::PotentialConfiguration& config)
        : Topology(std::move(particles)), config(config), graph_(std::move(graph)) {
    if (GraphTopology::graph().vertices().size() != GraphTopology::getNParticles()) {
        log::error("tried creating graph topology with {} vertices but only {} particles.",
                   GraphTopology::graph().vertices().size(), GraphTopology::getNParticles());
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

    std::unordered_map<api::BondType, std::vector<pot::BondConfiguration>, readdy::util::hash::EnumClassHash> bonds;
    std::unordered_map<api::AngleType, std::vector<pot::AngleConfiguration>, readdy::util::hash::EnumClassHash> angles;
    std::unordered_map<api::TorsionType, std::vector<pot::DihedralConfiguration>, readdy::util::hash::EnumClassHash> dihedrals;

    graph_.findNTuples([&](const graph_t::edge &tuple) {
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
            if (!v1->label().empty()) {
                ss << " (" << v1->label() << ")";
            }
            ss << " -- " << v2->particleIndex;
            if (!v2->label().empty()) {
                ss << " (" << v2->label() << ")";
            }
            ss << " has no bond configured! (See KernelContext.configureTopologyBondPotential())";

            throw std::invalid_argument(ss.str());
        }
    }, [&](const graph_t::path_len_2 &triple) {
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
    }, [&](const graph_t::path_len_3 &quadruple) {
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
                addBondedPotential(std::make_unique<harmonic_bond>(this, bond.second));
                break;
            };
        }
    }
    for (const auto &angle : angles) {
        switch (angle.first) {
            case api::AngleType::HARMONIC: {
                addAnglePotential(std::make_unique<harmonic_angle>(this, angle.second));
                break;
            };
        }
    }
    for (const auto &dih : dihedrals) {
        switch (dih.first) {
            case api::TorsionType::COS_DIHEDRAL: {
                addTorsionPotential(std::make_unique<cos_dihedral>(this, dih.second));
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

void GraphTopology::updateReactionRates() {
    _cumulativeRate = 0;
    for(auto&& reaction : reactions_) {
        const auto& r = std::get<0>(reaction);
        const auto rate = r.rate(*this);
        std::get<1>(reaction) = rate;
        _cumulativeRate += rate;
    }
    {
        std::stringstream ss;
        for(const auto& r : reactions_) {
            ss << std::get<1>(r) << " + ";
        }
        const void * address = static_cast<const void*>(this);
        std::stringstream ss2;
        ss2 << address;
        std::string name = ss2.str();
    }
}

void GraphTopology::addReaction(const reactions::TopologyReaction &reaction) {
    reactions_.push_back(std::make_tuple(reaction, 0));
}

void GraphTopology::addReaction(reactions::TopologyReaction &&reaction) {
    reactions_.push_back(std::make_tuple(std::move(reaction), 0));
}

const GraphTopology::topology_reactions &GraphTopology::registeredReactions() const {
    return reactions_;
}

GraphTopology::topology_reactions &GraphTopology::registeredReactions() {
    return reactions_;
}

std::vector<GraphTopology> GraphTopology::connectedComponents() {
    auto subGraphs = graph_.connectedComponentsDestructive();
    // generate particles list for each sub graph, update sub graph's vertices to obey this new list
    std::vector<particles_t> subGraphsParticles;
    {
        subGraphsParticles.reserve(subGraphs.size());
        for (auto itGraph = subGraphs.begin(); itGraph != subGraphs.end(); ++itGraph) {
            subGraphsParticles.emplace_back();
            auto &subParticles = subGraphsParticles.back();
            subParticles.reserve(itGraph->vertices().size());
            for (auto &vertex : itGraph->vertices()) {
                subParticles.emplace_back(particles.at(vertex.particleIndex));
                vertex.particleIndex = subParticles.size() - 1;
            }
        }
    }
    // create actual GraphTopology objects from graphs and particles
    std::vector<GraphTopology> components;
    {
        components.reserve(subGraphs.size());
        {
            auto it_graphs = subGraphs.begin();
            auto it_particles = subGraphsParticles.begin();
            for(; it_graphs != subGraphs.end(); ++it_graphs, ++it_particles) {
                components.emplace_back(std::move(*it_particles), std::move(*it_graphs), config);
                for(const auto& reaction : reactions_) {
                    components.back().addReaction(std::move(std::get<0>(reaction)));
                }
            }
        }
    }
    return std::move(components);
}

const bool GraphTopology::isDeactivated() const {
    return deactivated;
}

void GraphTopology::deactivate() {
    deactivated = true;
}

const GraphTopology::rate_t GraphTopology::cumulativeRate() const {
    return _cumulativeRate;
}

const bool GraphTopology::isNormalParticle(const Kernel &k) const {
    if(getNParticles() == 1){
        const auto particle_type = k.getKernelStateModel().getParticleType(particles.front());
        const auto& info = k.getKernelContext().particle_types().info_of(particle_type);
        return info.flavor != Particle::FLAVOR_TOPOLOGY;
    }
    return false;
}

}
}
}

