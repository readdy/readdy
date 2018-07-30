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
        if(!it_entries->is_deactivated() && it_entries->id == particle.id()) {
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




