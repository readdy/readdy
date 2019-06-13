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
 * @file SingleCPUParticleData.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.06.16
 */

#pragma once

#include <memory>
#include <vector>
#include <readdy/model/Particle.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace model {

//class SCPUParticleData;

struct Entry {
    using entries_vector = std::vector<Entry>;
    using particle_type = readdy::model::Particle;
    using force_type = particle_type::pos_type;
    using displacement_type = scalar;
    using topology_index_type = std::ptrdiff_t;

    explicit Entry(const particle_type &particle)
            : pos(particle.pos()), force(force_type()), type(particle.type()), deactivated(false),
              displacement(0), id(particle.id()) {}

    Entry(const Entry &) = default;

    Entry &operator=(const Entry &) = default;

    Entry(Entry &&) noexcept = default;

    Entry &operator=(Entry &&) noexcept = default;

    ~Entry() = default;

    bool is_deactivated() const {
        return deactivated;
    }

    const particle_type::pos_type &position() const {
        return pos;
    }

    force_type force;
    displacement_type displacement;
    particle_type::pos_type pos;
    topology_index_type topology_index {-1};
    particle_type::id_type id;
    particle_type::type_type type;
    bool deactivated;
};

template<typename EntryType>
class SCPUParticleData {
public:

    using entry_type = EntryType;
    using entries_vec = std::vector<EntryType>;
    using entry_index = typename entries_vec::size_type;
    using new_entries = std::vector<EntryType>;
    using particle_type = readdy::model::Particle;
    using top_particle_type = readdy::model::TopologyParticle;
    using force = particle_type::pos_type;
    using displacement = scalar;
    using iterator = typename entries_vec::iterator;
    using const_iterator = typename entries_vec::const_iterator;
    using entries_update = std::pair<new_entries, std::vector<entry_index>>;

    SCPUParticleData() = default;
    SCPUParticleData(const SCPUParticleData&) = delete;
    SCPUParticleData& operator=(const SCPUParticleData&) = delete;
    SCPUParticleData(SCPUParticleData&&) = default;
    SCPUParticleData& operator=(SCPUParticleData&&) = default;
    ~SCPUParticleData() = default;

    readdy::model::Particle getParticle(entry_index index) const {
        const auto& entry = *(entries.begin() + index);
        if(entry.deactivated) {
            log::error("Requested deactivated particle at index {}!", index);
        }
        return toParticle(entry);
    };

    readdy::model::Particle toParticle(const EntryType &e) const {
        return readdy::model::Particle(e.pos, e.type, e.id);
    }

    void addParticle(const particle_type &particle) {
        addParticles({particle});
    }

    void addParticles(const std::vector<particle_type> &particles) {
        for(const auto& p : particles) {
            if(!_blanks.empty()) {
                const auto idx = _blanks.back();
                _blanks.pop_back();
                entries.at(idx) = EntryType{p};
            } else {
                entries.emplace_back(p);
            }
        }
    }

    std::vector<typename entries_vec::size_type> addTopologyParticles(const std::vector<top_particle_type> &particles) {
        std::vector<typename entries_vec::size_type> indices;
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
    };

    void removeParticle(const particle_type &particle) {
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
    };

    void removeParticle(size_t index) {
        auto& p = *(entries.begin() + index);
        if(!p.deactivated) {
            _blanks.push_back(index);
            p.deactivated = true;
            // neighbors.at(index).clear();
        } else {
            log::error("Tried to remove particle (index={}), that was already removed!", index);
        }
    }

    iterator begin() {
        return entries.begin();
    }

    iterator end() {
        return entries.end();
    }

    const_iterator cbegin() const {
        return entries.cbegin();
    }

    const_iterator cend() const {
        return entries.cend();
    }

    const_iterator begin() const {
        return entries.begin();
    }

    const_iterator end() const {
        return entries.end();
    }

    EntryType &entry_at(entry_index idx) {
        return entries.at(idx);
    }

    entry_index size() const {
        return entries.size();
    }

    entry_index n_deactivated() const {
        return _blanks.size();
    }

    void reserve(std::size_t n) {
        entries.reserve(n);
    }

    void clear() {
        entries.clear();
        _blanks.clear();
    }

    const EntryType &entry_at(entry_index idx) const {
        return entries.at(idx);
    }

    const EntryType &centry_at(entry_index idx) const {
        return entries.at(idx);
    }

    entry_index addEntry(EntryType entry) {
        if(!_blanks.empty()) {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            entries.at(idx) = entry;
            return idx;
        }
        entries.push_back(entry);
        return entries.size()-1;
    }

    void removeEntry(entry_index idx) {
        auto &entry = entries.at(idx);
        if(!entry.is_deactivated()) {
            entry.deactivated = true;
            _blanks.push_back(idx);
        }
    };

    std::vector<entry_index> update(entries_update&& update_data) {
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
    };
    
    const std::vector<entry_index> &blanks() const {
        return _blanks;
    }

protected:

    entries_vec entries;
    std::vector<entry_index> _blanks;

};


}
}
}
}
