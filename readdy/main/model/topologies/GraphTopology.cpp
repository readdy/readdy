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

namespace readdy {
namespace model {
namespace top {

GraphTopology::GraphTopology(topology_type_type type,
                             const Topology::particle_indices &particles, const types_vec &types,
                             const model::Context& context, const model::StateModel *stateModel)
        : Topology(particles), _context(context), _topology_type(type), _stateModel(stateModel) {
    assert(types.size() == particles.size());
    std::size_t i = 0;
    for (auto itTypes = types.begin(); itTypes != types.end(); ++itTypes, ++i) {
        graph().addVertex(i, *itTypes);
    }
}

GraphTopology::GraphTopology(topology_type_type type,
                             Topology::particle_indices &&particles, graph::Graph &&graph,
                             const model::Context& context, const model::StateModel *stateModel)
        : Topology(std::move(particles)), _context(context), graph_(std::move(graph)), _topology_type(type),
          _stateModel(stateModel) {
    if (GraphTopology::graph().vertices().size() != GraphTopology::getNParticles()) {
        log::error("tried creating graph topology with {} vertices but only {} particles.",
                   GraphTopology::graph().vertices().size(), GraphTopology::getNParticles());
        throw std::invalid_argument("the number of particles and the number of vertices should match when creating"
                                            "a graph in this way!");
    }
    std::size_t idx = 0;
    auto &vertices = GraphTopology::graph().vertices();
    for (auto &vertex : vertices) {
        vertex.particleIndex = idx++;
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

    const auto &config = context().topology_registry().potential_configuration();

    graph_.findNTuples([&](const topology_graph::edge &tuple) {
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
            ss << " -- " << v2->particleIndex;
            ss << " has no bond configured! (See Context.configure_bond_potential())";

            throw std::invalid_argument(ss.str());
        }
    }, [&](const topology_graph::path_len_2 &triple) {
        const auto &v1 = std::get<0>(triple);
        const auto &v2 = std::get<1>(triple);
        const auto &v3 = std::get<2>(triple);
        auto it = config.anglePotentials.find(std::tie(v1->particleType(), v2->particleType(), v3->particleType()));
        if (it != config.anglePotentials.end()) {
            for (const auto &cfg : it->second) {
                angles[cfg.type].emplace_back(v2->particleIndex, v1->particleIndex, v3->particleIndex,
                                              cfg.forceConstant, cfg.equilibriumAngle);
            }
        }
    }, [&](const topology_graph::path_len_3 &quadruple) {
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
                addBondedPotential(std::make_unique<harmonic_bond>(bond.second));
                break;
            };
        }
    }
    for (const auto &angle : angles) {
        switch (angle.first) {
            case api::AngleType::HARMONIC: {
                addAnglePotential(std::make_unique<harmonic_angle>(angle.second));
                break;
            };
        }
    }
    for (const auto &dih : dihedrals) {
        switch (dih.first) {
            case api::TorsionType::COS_DIHEDRAL: {
                addTorsionPotential(std::make_unique<cos_dihedral>(dih.second));
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

void GraphTopology::updateReactionRates(const TopologyRegistry::structural_reactions &reactions) {
    _cumulativeRate = 0;
    _reaction_rates.resize(reactions.size());
    auto it = _reaction_rates.begin();
    for(auto&& reaction : reactions) {
        const auto rate = reaction.rate(*this);
        *it = rate;
        _cumulativeRate += rate;
        ++it;
    }
}

std::vector<GraphTopology> GraphTopology::connectedComponents() {
    auto subGraphs = graph_.connectedComponentsDestructive();
    // generate particles list for each sub graph, update sub graph's vertices to obey this new list
    std::vector<particle_indices> subGraphsParticles;
    {
        subGraphsParticles.reserve(subGraphs.size());
        for (auto &subGraph : subGraphs) {
            subGraphsParticles.emplace_back();
            auto &subParticles = subGraphsParticles.back();
            subParticles.reserve(subGraph.vertices().size());
            for (auto &vertex : subGraph.vertices()) {
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
                components.emplace_back(_topology_type, std::move(*it_particles), std::move(*it_graphs), _context, _stateModel);
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

const GraphTopology::topology_reaction_rate GraphTopology::cumulativeRate() const {
    return _cumulativeRate;
}

const bool GraphTopology::isNormalParticle(const Kernel &k) const {
    if(getNParticles() == 1){
        const auto particle_type = k.stateModel().getParticleType(particles.front());
        const auto& info = k.context().particle_types().info_of(particle_type);
        return info.flavor != particleflavor::NORMAL;
    }
    return false;
}

void GraphTopology::appendParticle(particle_index newParticle, particle_type_type newParticleType,
                                   particle_index counterPart, particle_type_type counterPartType) {
    auto it = std::find(particles.begin(), particles.end(), counterPart);
    if(it != particles.end()) {
        auto counterPartIdx = std::distance(particles.begin(), it);

        particles.push_back(newParticle);
        graph().addVertex(particles.size() - 1, newParticleType);

        auto newParticleIt = std::prev(graph().vertices().end());
        auto otherParticleIt = std::next(graph().vertices().begin(), counterPartIdx);
        otherParticleIt->setParticleType(counterPartType);

        graph().addEdge(newParticleIt, otherParticleIt);
    } else {
        log::critical("counterPart {} was not contained in topology, this should not happen", counterPart);
    }
}

topology_type_type &GraphTopology::type() {
    return _topology_type;
}

const topology_type_type &GraphTopology::type() const {
    return _topology_type;
}

const GraphTopology::topology_reaction_rates &GraphTopology::rates() const {
    return _reaction_rates;
}

void GraphTopology::appendTopology(GraphTopology &other, Topology::particle_index otherParticle,
                                   particle_type_type otherNewParticleType, Topology::particle_index thisParticle,
                                   particle_type_type thisNewParticleType, topology_type_type newType) {
    auto &otherGraph = other.graph();
    auto &thisGraph = graph();

    auto former_begin = otherGraph.vertices().begin();
    auto former_n_vertices = particles.size();

    auto other_vert = other.vertexForParticle(otherParticle);
    auto this_vert = vertexForParticle(thisParticle);

    // insert other particles into this' particles
    particles.insert(std::end(particles), std::begin(other.particles), std::end(other.particles));
    // move other graph into this graph
    thisGraph.vertices().splice(thisGraph.vertices().end(), otherGraph.vertices());

    for(auto it = former_begin; it != graph().vertices().end(); ++it) {
        it->particleIndex = former_n_vertices;
        ++former_n_vertices;
    }

    // add edge between the formerly two topologies
    graph().addEdge(other_vert, this_vert);
    other_vert->setParticleType(otherNewParticleType);
    this_vert->setParticleType(thisNewParticleType);

    _topology_type = newType;
}

graph::Graph::vertex_ref GraphTopology::vertexForParticle(Topology::particle_index particle) {
    auto it = std::find(particles.begin(), particles.end(), particle);
    if(it != particles.end()) {
        return std::next(graph_.vertices().begin(), std::distance(particles.begin(), it));
    }
    return graph_.vertices().end();
}

const model::Context &GraphTopology::context() const {
    return _context.get();
}

std::vector<Particle> GraphTopology::fetchParticles() const {
    if(!_stateModel) {
        throw std::logic_error("Cannot fetch particles if no state model was provided!");
    }
    return _stateModel->getParticlesForTopology(*this);
}

Particle GraphTopology::particleForVertex(graph::Graph::vertex_ref vertexRef) const {
    return particleForVertex(*vertexRef);
}

Particle GraphTopology::particleForVertex(const vertex &vertex) const {
    if(!_stateModel) {
        throw std::logic_error("Cannot fetch particle if not state model was provided!");
    }
    return _stateModel->getParticleForIndex(vertex.particleIndex);
}


}
}
}
