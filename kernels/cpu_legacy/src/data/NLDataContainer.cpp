/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file NLDataContainer.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 14.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/cpu/data/NLDataContainer.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace data {

void NLDataContainer::reserve(std::size_t n) {
    _entries.reserve(n);
    _neighbors.reserve(n);
    _displacements.reserve(n);
}

NLDataContainer::size_type NLDataContainer::addEntry(Entry &&entry) {
    if (!_blanks.empty()) {
        const auto idx = _blanks.back();
        _blanks.pop_back();
        _entries.at(idx) = std::move(entry);
        _neighbors.at(idx).clear();
        _displacements.at(idx) = 0;
        return idx;
    }

    _entries.push_back(std::move(entry));
    _neighbors.emplace_back();
    _displacements.emplace_back();
    return _entries.size() - 1;
}

void NLDataContainer::addParticles(const std::vector<Particle> &particles) {
    for (const auto &p : particles) {
        if (!_blanks.empty()) {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _entries.at(idx) = Entry(p);
            _neighbors.at(idx).clear();
            _displacements.at(idx) = 0;
        } else {
            _entries.emplace_back(p);
            _neighbors.emplace_back();
            _displacements.emplace_back();
        }
    }
}

std::vector<NLDataContainer::size_type>
NLDataContainer::addTopologyParticles(const std::vector<TopologyParticle> &topologyParticles) {
    std::vector<size_type> indices;
    indices.reserve(topologyParticles.size());
    for (const auto &p : topologyParticles) {
        if (!_blanks.empty()) {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _entries.at(idx) = Entry(p);
            _neighbors.at(idx).clear();
            _displacements.at(idx) = 0;
            indices.push_back(idx);
        } else {
            _entries.emplace_back(p);
            _neighbors.emplace_back();
            _displacements.emplace_back();
            indices.push_back(_entries.size() - 1);
        }
    }
    return indices;
}

std::vector<NLDataContainer::size_type> NLDataContainer::update(DataUpdate &&update) {
    std::vector<size_type> result;

    auto &&newEntries = std::move(std::get<0>(update));
    auto &&removedEntries = std::move(std::get<1>(update));
    result.reserve(newEntries.size());

    auto it_del = removedEntries.begin();
    for (auto &&newEntry : newEntries) {
        if (it_del != removedEntries.end()) {
            _entries.at(*it_del) = std::move(newEntry);
            _neighbors.at(*it_del).clear();
            _displacements.at(*it_del) = 0;
            result.push_back(*it_del);
            ++it_del;
        } else {
            result.push_back(addEntry(std::move(newEntry)));
        }
    }
    while (it_del != removedEntries.end()) {
        removeEntry(*it_del);
        ++it_del;
    }

    return result;
}

void NLDataContainer::displace(size_type ix, const Particle::pos_type &delta) {
    auto &entry = _entries.at(ix);
    entry.pos += delta;
    _context.get().fixPositionFun()(entry.pos);
    _displacements.at(ix) += std::sqrt(delta * delta);
}

NLDataContainer::Neighbors &NLDataContainer::neighbors_at(DataContainer::size_type index) {
    return _neighbors.at(index);
}

const NLDataContainer::Neighbors &NLDataContainer::neighbors_at(DataContainer::size_type index) const {
    return _neighbors.at(index);
}

const NLDataContainer::Neighbors &NLDataContainer::cneighbors_at(DataContainer::size_type index) const {
    return _neighbors.at(index);
}

const NLDataContainer::NeighborList &NLDataContainer::neighbors() const {
    return _neighbors;
}

NLDataContainer::NeighborList &NLDataContainer::neighbors() {
    return _neighbors;
}

NLDataContainer::NLDataContainer(const model::Context &context, const util::thread::Config &threadConfig)
        : DataContainer(context, threadConfig), _neighbors{} {
}

NLDataContainer::NLDataContainer(EntryDataContainer *data) : NLDataContainer(data->context(), data->threadConfig()) {
    _entries = data->entries();
    _blanks = data->blanks();
    _neighbors.resize(_entries.size());
    _displacements.resize(_entries.size());
}

std::vector<scalar> &NLDataContainer::displacements() {
    return _displacements;
}

const std::vector<scalar> &NLDataContainer::displacements() const {
    return _displacements;
}

}
}
}
}
