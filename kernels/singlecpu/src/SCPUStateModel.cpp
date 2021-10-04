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
 * @file SCPUStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include <algorithm>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>

namespace readdy::kernel::scpu {

SCPUStateModel::SCPUStateModel(const readdy::model::Context &context, topology_action_factory const *const taf)
        : _context(context) {
    neighborList = std::make_unique<model::CellLinkedList>(particleData, context);
    topologyActionFactory = taf;
}

std::vector<Vec3> readdy::kernel::scpu::SCPUStateModel::getParticlePositions() const {
    const auto &data = particleData;
    std::vector<Vec3> target{};
    target.reserve(data.size());
    for (const auto &entry : data) {
        if (!entry.is_deactivated()) target.push_back(entry.position());
    }
    return target;
}

std::vector<readdy::model::Particle> SCPUStateModel::getParticles() const {
    const auto &data = particleData;
    std::vector<readdy::model::Particle> result;
    result.reserve(data.size());
    for (const auto &entry : data) {
        if (!entry.is_deactivated()) {
            result.push_back(data.toParticle(entry));
        }
    }
    return result;
}


readdy::model::top::GraphTopology *const SCPUStateModel::addTopology(TopologyTypeId type, const std::vector<readdy::model::Particle> &particles) {
    std::vector<std::size_t> indices = particleData.addTopologyParticles(particles);
    readdy::model::top::Graph graph;
    for(auto index : indices) {
        graph.addVertex(readdy::model::top::VertexData{index});
    }

    auto it = _topologies.emplace_back(
            std::make_unique<topology>(type, std::move(graph), _context.get(), this)
    );
    const auto idx = std::distance(topologies().begin(), it);
    for(auto v : (*it)->graph().vertices()) {
        if(!v.deactivated()) {
            particleData.entry_at(v->particleIndex).topology_index = idx;
        }
    }
    return it->get();
}

std::vector<readdy::model::top::GraphTopology*> SCPUStateModel::getTopologies() {
    std::vector<readdy::model::top::GraphTopology*> result;
    result.reserve(_topologies.size() - _topologies.n_deactivated());
    for(const auto& top : _topologies) {
        if(!top->isDeactivated()) {
            result.push_back(top.get());
        }
    }
    return result;
}

void SCPUStateModel::insert_topology(SCPUStateModel::topology &&top) {
    auto it = _topologies.push_back(std::make_unique<topology>(std::move(top)));
    auto idx = std::distance(_topologies.begin(), it);
    const auto& vertices = it->get()->graph().vertices();
    auto& data = particleData;
    std::for_each(vertices.begin(), vertices.end(), [idx, &data](const auto& vertex) {
        data.entry_at(vertex->particleIndex).topology_index = idx;
    });
}

void SCPUStateModel::resetReactionCounts() {
    if(!reactionCounts().empty()) {
        for(auto &e : reactionCounts()) {
            e.second = 0;
        }
    } else {
        const auto &reactions = _context.get().reactions();
        for (const auto &entry : reactions.order1()) {
            for (auto reaction : entry.second) {
                reactionCounts()[reaction->id()] = 0;
            }
        }
        for (const auto &entry : reactions.order2()) {
            for (auto reaction : entry.second) {
                reactionCounts()[reaction->id()] = 0;
            }
        }
    }
}

void SCPUStateModel::toDenseParticleIndices(std::vector<std::size_t>::iterator begin,
                                            std::vector<std::size_t>::iterator end) const {
    const auto &blanks = particleData.blanks();
    std::transform(begin, end, begin, [&blanks](const std::size_t &ix) {
        auto result = ix;
        for(auto blankIx : blanks) {
            if(blankIx < ix) {
                --result;
            }
        }
        return result;
    });
}

void SCPUStateModel::clear() {
    particleData.clear();
    _topologies.clear();
    reactionRecords().clear();
    resetReactionCounts();
    virial() = {};
    energy() = 0;
}

}
