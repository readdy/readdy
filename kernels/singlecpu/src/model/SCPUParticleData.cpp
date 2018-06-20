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

#include <memory>

#include <readdy/common/logging.h>
#include <readdy/kernel/singlecpu/model/SCPUParticleData.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace model {

readdy::model::Particle SCPUParticleData::getParticle(const entry_index index) const {
    const auto& entry = *(entries.begin() + index);
    if(entry.deactivated) {
        log::error("Requested deactivated particle at index {}!", index);
    }
    return toParticle(entry);
}

std::vector<SCPUParticleData::entries_vec::size_type>
SCPUParticleData::addTopologyParticles(const std::vector<SCPUParticleData::top_particle_type> &particles) {
    std::vector<entries_vec::size_type> indices;
    indices.reserve(particles.size());
    for(const auto& p : particles) {
        if(!_blanks.empty()) {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            entries.at(idx) = Entry{p};
            indices.push_back(idx);
        } else {
            indices.push_back(entries.size());
            entries.emplace_back(p);
        }
    }
    return indices;
}

void SCPUParticleData::removeParticle(const SCPUParticleData::particle_type &particle) {
    auto it_entries = begin();
    std::size_t idx = 0;
    for(; it_entries != end(); ++it_entries, ++idx) {
        if(!it_entries->is_deactivated() && it_entries->id == particle.getId()) {
            _blanks.push_back(idx);
            it_entries->deactivated = true;
            return;
        }
    }
    log::error("Tried to remove particle ({}) which did not exist or was already deactivated!", particle);
}

std::vector<SCPUParticleData::entry_index> SCPUParticleData::update(SCPUParticleData::entries_update &&update_data) {
    std::vector<entry_index> result;

    auto &&newEntries = std::move(std::get<0>(update_data));
    auto &&removedEntries = std::move(std::get<1>(update_data));
    result.reserve(newEntries.size());

    auto it_del = removedEntries.begin();
    for(const auto& newEntry : newEntries) {
        if(it_del != removedEntries.end()) {
            entries.at(*it_del) = newEntry;
            result.push_back(*it_del);
            ++it_del;
        } else {
            result.push_back(addEntry(newEntry));
        }
    }
    while(it_del != removedEntries.end()) {
        removeEntry(*it_del);
        ++it_del;
    }

    return result;
}

void SCPUParticleData::removeEntry(SCPUParticleData::entry_index idx) {
    auto &entry = entries.at(idx);
    if(!entry.is_deactivated()) {
        entry.deactivated = true;
        _blanks.push_back(idx);
    }
}

}
}
}
}




