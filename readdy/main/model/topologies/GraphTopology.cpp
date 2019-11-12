/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file GraphTopology.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.03.17
 * @copyright BSD-3
 */

#include <sstream>

#include <readdy/model/Kernel.h>
#include <readdy/model/topologies/GraphTopology.h>


namespace readdy::model::top {

GraphTopology::GraphTopology(TopologyTypeId type, Graph graph,
                             const model::Context& context, const model::StateModel *stateModel)
        : Topology(), _context(context), _topology_type(type), _stateModel(stateModel), _cumulativeRate(0),
        _graph(std::move(graph)) {}

void GraphTopology::configure() {
    validate();

    bondedPotentials.clear();
    anglePotentials.clear();
    torsionPotentials.clear();

    std::unordered_map<api::BondType, std::vector<pot::BondConfiguration>, readdy::util::hash::EnumClassHash> bonds;
    std::unordered_map<api::AngleType, std::vector<pot::AngleConfiguration>, readdy::util::hash::EnumClassHash> angles;
    std::unordered_map<api::TorsionType, std::vector<pot::DihedralConfiguration>, readdy::util::hash::EnumClassHash> dihedrals;

    const auto &config = context().topologyRegistry().potentialConfiguration();

    _graph.findNTuples([&](const Graph::Edge &tuple) {
        auto [i1, i2] = tuple;
        const auto& v1 = _graph.vertices().at(i1);
        const auto& v2 = _graph.vertices().at(i2);
        auto it = config.pairPotentials.find(std::make_tuple(typeOf(v1), typeOf(v2)));
        if (it != config.pairPotentials.end()) {
            for (const auto &cfg : it->second) {
                bonds[cfg.type].emplace_back(v1->particleIndex, v2->particleIndex, cfg.forceConstant, cfg.length);
            }
        } else {
            std::ostringstream ss;
            auto p1 = particleForVertex(v1);
            auto p2 = particleForVertex(v2);

            ss << "The edge " << v1->particleIndex << " (" << context().particleTypes().nameOf(p1.type()) << ")";
            ss << " -- " << v2->particleIndex << " (" << context().particleTypes().nameOf(p2.type()) << ")";
            ss << " has no bond configured!";

            throw std::invalid_argument(ss.str());
        }
    }, [&](const Graph::Path3 &triple) {
        auto [i1, i2, i3] = triple;
        const auto& v1 = _graph.vertices().at(i1);
        const auto& v2 = _graph.vertices().at(i2);
        const auto& v3 = _graph.vertices().at(i3);
        auto it = config.anglePotentials.find(std::make_tuple(typeOf(v1), typeOf(v2), typeOf(v3)));
        if (it != config.anglePotentials.end()) {
            for (const auto &cfg : it->second) {
                angles[cfg.type].emplace_back(v1->particleIndex, v2->particleIndex, v3->particleIndex,
                                              cfg.forceConstant, cfg.equilibriumAngle);
            }
        }
    }, [&](const Graph::Path4 &quadruple) {
        auto [i1, i2, i3, i4] = quadruple;
        const auto& v1 = _graph.vertices().at(i1);
        const auto& v2 = _graph.vertices().at(i2);
        const auto& v3 = _graph.vertices().at(i3);
        const auto& v4 = _graph.vertices().at(i4);
        auto it = config.torsionPotentials.find(std::make_tuple(typeOf(v1), typeOf(v2), typeOf(v3), typeOf(v4)));
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
                addBondedPotential(std::make_unique<HarmonicBond>(bond.second));
                break;
            }
        }
    }
    for (const auto &angle : angles) {
        switch (angle.first) {
            case api::AngleType::HARMONIC: {
                addAnglePotential(std::make_unique<HarmonicAngle>(angle.second));
                break;
            }
        }
    }
    for (const auto &dih : dihedrals) {
        switch (dih.first) {
            case api::TorsionType::COS_DIHEDRAL: {
                addTorsionPotential(std::make_unique<CosineDihedral>(dih.second));
                break;
            }
        }
    }
}

std::vector<GraphTopology> GraphTopology::connectedComponents() {
    auto subGraphs = _graph.connectedComponents();
    // create actual GraphTopology objects from graphs and particles
    std::vector<GraphTopology> components;
    {
        components.reserve(subGraphs.size());
        {
            auto it_graphs = subGraphs.begin();
            for(; it_graphs != subGraphs.end(); ++it_graphs) {
                components.emplace_back(_topology_type, std::move(*it_graphs), _context, _stateModel);
            }
        }
    }
    return std::move(components);
}

bool GraphTopology::isNormalParticle(const Kernel &k) const {
    if(_graph.nVertices() == 1){
        for(const auto &v : _graph.vertices()) {
            if (!v.deactivated()) {
                const auto particle_type = k.stateModel().getParticleType(typeOf(v));
                const auto& info = k.context().particleTypes().infoOf(particle_type);
                return info.flavor == particleflavor::NORMAL;
            }
        }
        throw std::logic_error("graph has 1 vertex but only tombstone entries");
    }
    return false;
}

void GraphTopology::appendParticle(VertexData::ParticleIndex newParticle, ParticleTypeId newParticleType,
                                   VertexData::ParticleIndex counterPart, ParticleTypeId counterPartType) {
    auto it = vertexIndexForParticle(counterPart);
    if(it == _graph.vertices().end()) {
        throw std::invalid_argument(fmt::format("Vertex with particle {} was not in graph or deactivated!", counterPart));
    }
    auto ix = _graph.addVertex(VertexData{
        .particleIndex = newParticle,
    });

    (*it)->particleType = counterPartType;
    _graph.addEdge(ix, std::distance(_graph.vertices().begin(), it));
}

void GraphTopology::appendTopology(GraphTopology &other, VertexData::ParticleIndex otherParticle,
                                   ParticleTypeId otherNewParticleType, VertexData::ParticleIndex thisParticle,
                                   ParticleTypeId thisNewParticleType, TopologyTypeId newType) {
    auto &otherGraph = other.graph();

    if(!otherGraph.vertices().empty()) {

        auto &thisGraph = graph();

        auto it = vertexIndexForParticle(thisParticle);
        if(it == _graph.vertices().end()) {
            throw std::invalid_argument("...");
        }
        (*it)->particleType = thisNewParticleType;

        auto itOther = other.vertexIndexForParticle(otherParticle);
        if(itOther == other.graph().vertices().end()) {
            throw std::invalid_argument("...");
        }
        (*itOther)->particleType = otherNewParticleType;

        thisGraph.append(otherGraph, std::distance(_graph.vertices().begin(), it),
                         std::distance(otherGraph.vertices().begin(), itOther));
        _topology_type = newType;
    } else {
        log::warn("encountered empty topology which was deactivated={}", other.isDeactivated());
    }
}

std::vector<Particle> GraphTopology::fetchParticles() const {
    if(!_stateModel) {
        throw std::logic_error("Cannot fetch particles if state model was not provided!");
    }
    return _stateModel->getParticlesForTopology(*this);
}

typename Graph::VertexList::iterator GraphTopology::vertexIndexForParticle(VertexData::ParticleIndex index) {
    return std::find_if(_graph.vertices().begin(), _graph.vertices().end(), [index](const auto& v) {
        return !v.deactivated() && v->particleIndex == index;
    });
}

Particle GraphTopology::particleForVertex(const Vertex &vertex) const {
    if(!_stateModel) {
        throw std::logic_error("Cannot fetch particle if state model was not provided!");
    }
    return _stateModel->getParticleForIndex(vertex->particleIndex);
}

ParticleTypeId GraphTopology::typeOf(Graph::VertexIndex vertex) const {
    return typeOf(graph().vertices().at(vertex));
}

ParticleTypeId GraphTopology::typeOf(const Vertex &v) const {
    return particleForVertex(v).type();
}

}
