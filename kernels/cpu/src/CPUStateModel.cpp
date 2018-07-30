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
 * @file CPUStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12/11/17
 */


#include <future>
#include <readdy/kernel/cpu/CPUStateModel.h>

namespace readdy {
namespace kernel {
namespace cpu {

namespace thd = readdy::util::thread;

using entries_it = CPUStateModel::data_type::Entries::iterator;
using topologies_it = std::vector<std::unique_ptr<readdy::model::top::GraphTopology>>::const_iterator;
using pot1Map = readdy::model::potentials::PotentialRegistry::PotentialsO1Map;
using pot2Map = readdy::model::potentials::PotentialRegistry::PotentialsO2Map;
using dist_fun = readdy::model::Context::shortest_dist_fun;

const std::vector<Vec3> CPUStateModel::getParticlePositions() const {
    const auto data = getParticleData();
    std::vector<Vec3> target{};
    target.reserve(data->size());
    for (const auto &entry : *data) {
        if (!entry.deactivated) target.push_back(entry.pos);
    }
    return target;
}

const std::vector<readdy::model::Particle> CPUStateModel::getParticles() const {
    const auto data = getParticleData();
    std::vector<readdy::model::Particle> result;
    result.reserve(data->size());
    for (const auto &entry : *data) {
        if (!entry.deactivated) {
            result.push_back(data->toParticle(entry));
        }
    }
    return result;
}

CPUStateModel::CPUStateModel(data_type &data, const readdy::model::Context &context,
                             thread_pool &pool,
                             readdy::model::top::TopologyActionFactory const *const taf)
        : _pool(pool), _context(context), _topologyActionFactory(*taf), _data(data) {
    _neighborList = std::make_unique<neighbor_list>(_data.get(), _context.get(), _pool.get());
    _reorderConnection = std::make_unique<readdy::signals::scoped_connection>(
            getParticleData()->registerReorderEventListener([this](const std::vector<std::size_t> &indices) -> void {
                for (auto &top : _topologies) {
                    if(!top->isDeactivated()) top->permuteIndices(indices);
                }
            }));
}

readdy::model::top::GraphTopology *const
CPUStateModel::addTopology(TopologyTypeId type, const std::vector<readdy::model::TopologyParticle> &particles) {
    std::vector<std::size_t> ids = getParticleData()->addTopologyParticles(particles);
    std::vector<ParticleTypeId> types;
    types.reserve(ids.size());
    for (const auto &p : particles) {
        types.push_back(p.type());
    }
    auto it = _topologies.push_back(std::make_unique<topology>(type, std::move(ids), std::move(types), _context.get(), this));
    const auto idx = std::distance(topologies().begin(), it);
    for(const auto p : (*it)->getParticles()) {
        getParticleData()->entry_at(p).topology_index = idx;
    }
    return it->get();
}

std::vector<readdy::model::top::GraphTopology*> CPUStateModel::getTopologies() {
    std::vector<readdy::model::top::GraphTopology*> result;
    result.reserve(_topologies.size() - _topologies.n_deactivated());
    for(const auto& top : _topologies) {
        if(!top->isDeactivated()) {
            result.push_back(top.get());
        }
    }
    return result;
}

const readdy::model::top::GraphTopology *CPUStateModel::getTopologyForParticle(readdy::model::top::Topology::particle_index particle) const {
    const auto& entry = getParticleData()->entry_at(particle);
    if(!entry.deactivated) {
        if(entry.topology_index >= 0) {
            return _topologies.at(static_cast<topologies_vec::size_type>(entry.topology_index)).get();
        }
        log::trace("requested particle {} of type {} had no assigned topology", particle, entry.type);
        return nullptr;
    }
    throw std::logic_error(fmt::format("requested particle was deactivated in getTopologyForParticle(p={})", particle));
}

readdy::model::top::GraphTopology *CPUStateModel::getTopologyForParticle(readdy::model::top::Topology::particle_index particle) {
    const auto& entry = getParticleData()->entry_at(particle);
    if(!entry.deactivated) {
        if(entry.topology_index >= 0) {
            return _topologies.at(static_cast<topologies_vec::size_type>(entry.topology_index)).get();
        }
        log::trace("requested particle {} of type {} had no assigned topology", particle, entry.type);
        return nullptr;
    }
    throw std::logic_error(fmt::format("requested particle was deactivated in getTopologyForParticle(p={})", particle));
}

void CPUStateModel::insert_topology(CPUStateModel::topology &&top) {
    auto it = _topologies.push_back(std::make_unique<topology>(std::move(top)));
    auto idx = std::distance(_topologies.begin(), it);
    const auto& particles = it->get()->getParticles();
    auto& data = *getParticleData();
    std::for_each(particles.begin(), particles.end(), [idx, &data](const topology::particle_index p) {
        data.entry_at(p).topology_index = idx;
    });
}

void CPUStateModel::resetReactionCounts() {
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

void CPUStateModel::toDenseParticleIndices(std::vector<std::size_t>::iterator begin,
                                           std::vector<std::size_t>::iterator end) const {
    const auto &blanks = _data.get().blanks();
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

void CPUStateModel::clear() {
    getParticleData()->clear();
    topologies().clear();
    reactionRecords().clear();
    resetReactionCounts();
    virial() = {};
    energy() = 0;
}


}
}
}