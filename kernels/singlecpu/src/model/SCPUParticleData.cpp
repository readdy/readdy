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
 * @file SingleCPUParticleData.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#include <readdy/common/make_unique.h>
#include <readdy/common/logging.h>
#include <readdy/kernel/singlecpu/model/SCPUParticleData.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace model {

SCPUParticleData::iterator SCPUParticleData::begin() {
    return entries.begin();
}

SCPUParticleData::iterator SCPUParticleData::end() {
    return entries.end();
}

SCPUParticleData::const_iterator SCPUParticleData::cbegin() const {
    return entries.cbegin();
}

SCPUParticleData::const_iterator SCPUParticleData::cend() const {
    return entries.cend();
}

SCPUParticleData::const_iterator SCPUParticleData::begin() const {
    return entries.begin();
}

SCPUParticleData::const_iterator SCPUParticleData::end() const {
    return entries.end();
}

readdy::model::Particle SCPUParticleData::getParticle(const index_t index) const {
    const auto& entry = *(entries.begin() + index);
    if(entry.deactivated) {
        log::console()->error("Requested deactivated particle at index {}!", index);
    }
    return toParticle(entry);
}

readdy::model::Particle SCPUParticleData::toParticle(const Entry &e) const {
    return readdy::model::Particle(e.pos, e.type, e.id);
}

void SCPUParticleData::addParticle(const SCPUParticleData::particle_type &particle) {
    addParticles({particle});
}

void SCPUParticleData::addParticles(const std::vector<SCPUParticleData::particle_type> &particles) {
    for(const auto& p : particles) {
        if(!blanks.empty()) {
            const auto idx = blanks.back();
            blanks.pop_back();
            entries.at(idx) = {p};
        } else {
            entries.push_back({p});
        }
    }
}

std::vector<SCPUParticleData::entries_t::size_type>
SCPUParticleData::addTopologyParticles(const std::vector<SCPUParticleData::top_particle_type> &particles) {
    std::vector<entries_t::size_type> indices;
    indices.reserve(particles.size());
    for(const auto& p : particles) {
        if(!blanks.empty()) {
            const auto idx = blanks.back();
            blanks.pop_back();
            entries.at(idx) = {p};
            indices.push_back(idx);
        } else {
            indices.push_back(entries.size());
            entries.push_back({p});
        }
    }
    return indices;
}

void SCPUParticleData::removeParticle(const SCPUParticleData::particle_type &particle) {
    auto it_entries = begin();
    std::size_t idx = 0;
    for(; it_entries != end(); ++it_entries, ++idx) {
        if(!it_entries->is_deactivated() && it_entries->id == particle.getId()) {
            blanks.push_back(idx);
            it_entries->deactivated = true;
            return;
        }
    }
    log::console()->error("Tried to remove particle ({}) which did not exist or was already deactivated!", particle);
}

void SCPUParticleData::removeParticle(const size_t index) {
    auto& p = *(entries.begin() + index);
    if(!p.deactivated) {
        blanks.push_back(index);
        p.deactivated = true;
        // neighbors.at(index).clear();
    } else {
        log::console()->error("Tried to remove particle (index={}), that was already removed!", index);
    }
}

Entry &SCPUParticleData::entry_at(SCPUParticleData::index_t idx) {
    return entries.at(idx);
}

const Entry &SCPUParticleData::entry_at(SCPUParticleData::index_t idx) const {
    return entries.at(idx);
}

const Entry &SCPUParticleData::centry_at(SCPUParticleData::index_t idx) const {
    return entries.at(idx);
}

SCPUParticleData::index_t SCPUParticleData::size() const {
    return entries.size();
}

void SCPUParticleData::clear() {
    entries.clear();
    blanks.clear();
}

std::vector<SCPUParticleData::index_t> SCPUParticleData::update(SCPUParticleData::update_t &&update_data) {
    std::vector<index_t> result;

    auto &&newEntries = std::move(std::get<0>(update_data));
    auto &&removedEntries = std::move(std::get<1>(update_data));
    result.reserve(newEntries.size());

    auto it_del = removedEntries.begin();
    for(auto&& newEntry : newEntries) {
        if(it_del != removedEntries.end()) {
            entries.at(*it_del) = std::move(newEntry);
            result.push_back(*it_del);
            ++it_del;
        } else {
            result.push_back(addEntry(std::move(newEntry)));
        }
    }
    while(it_del != removedEntries.end()) {
        removeEntry(*it_del);
        ++it_del;
    }

    return result;
}

SCPUParticleData::index_t SCPUParticleData::addEntry(Entry &&entry) {
    if(!blanks.empty()) {
        const auto idx = blanks.back();
        blanks.pop_back();
        entries.at(idx) = std::move(entry);
        return idx;
    } else {
        entries.push_back(std::move(entry));
        return entries.size()-1;
    }
}

void SCPUParticleData::removeEntry(SCPUParticleData::index_t idx) {
    auto &entry = entries.at(idx);
    if(!entry.is_deactivated()) {
        entry.deactivated = true;
        blanks.push_back(idx);
    }
}

SCPUParticleData::index_t SCPUParticleData::n_deactivated() const {
    return blanks.size();
}

bool Entry::is_deactivated() const {
    return deactivated;
}

const readdy::model::Particle::pos_type &Entry::position() const {
    return pos;
}


}
}
}
}




