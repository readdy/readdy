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
 * @file DefaultDataContainer.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 14.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/cpu_legacy/data/DefaultDataContainer.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace data {

DefaultDataContainer::DefaultDataContainer(const readdy::model::Context &context, const util::thread::Config &threadConfig)
        : DataContainer(context, threadConfig) {}

DefaultDataContainer::DefaultDataContainer(EntryDataContainer *entryDataContainer)
        : DataContainer(entryDataContainer->context(), entryDataContainer->threadConfig()) {
    _entries = entryDataContainer->entries();
    //_entries.insert(_entries.end(), entryDataContainer.entries().begin(), entryDataContainer.entries().end());
    _blanks = entryDataContainer->blanks();
    // _blanks.insert(_blanks.end(), entryDataContainer.blanks().begin(), entryDataContainer.blanks().end());
}

void DefaultDataContainer::reserve(std::size_t n) {
    _entries.reserve(n);
}

DefaultDataContainer::size_type DefaultDataContainer::addEntry(Entry &&entry) {
    if(!_blanks.empty()) {
        const auto idx = _blanks.back();
        _blanks.pop_back();
        _entries.at(idx) = std::move(entry);
        return idx;
    }

    _entries.push_back(std::move(entry));
    return _entries.size()-1;
}

void DefaultDataContainer::addParticles(const std::vector<Particle> &particles) {
    for(const auto& p : particles) {
        if(!_blanks.empty()) {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _entries.at(idx) = Entry(p);
        } else {
            _entries.emplace_back(p);
        }
    }
}

std::vector<DefaultDataContainer::size_type>
DefaultDataContainer::addTopologyParticles(const std::vector<TopologyParticle> &topologyParticles) {
    std::vector<size_type> indices;
    indices.reserve(topologyParticles.size());
    for(const auto& p : topologyParticles) {
        if(!_blanks.empty()) {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _entries.at(idx) = Entry(p);
            indices.push_back(idx);
        } else {
            _entries.emplace_back(p);
            indices.push_back(_entries.size()-1);
        }
    }
    return indices;
}

std::vector<DefaultDataContainer::size_type> DefaultDataContainer::update(DataUpdate &&update) {
    auto &&newEntries = std::move(std::get<0>(update));
    auto &&removedEntries = std::move(std::get<1>(update));

    auto it_del = removedEntries.begin();
    for(auto&& newEntry : newEntries) {
        if(it_del != removedEntries.end()) {
            _entries.at(*it_del) = std::move(newEntry);
            ++it_del;
        } else {
            addEntry(std::move(newEntry));
        }
    }
    while(it_del != removedEntries.end()) {
        removeEntry(*it_del);
        ++it_del;
    }
    return {};
}

void DefaultDataContainer::displace(size_type index, const Particle::pos_type &delta) {
    auto &entry = _entries.at(index);
    entry.pos += delta;
    _context.get().fixPositionFun()(entry.pos);
}



}
}
}
}
