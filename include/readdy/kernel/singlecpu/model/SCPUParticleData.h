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

class SCPUParticleData;

struct Entry {
    using entries_vector = std::vector<Entry>;
    using particle_type = readdy::model::Particle;
    using force_type = particle_type::pos_type;
    using displacement_type = scalar;
    using topology_index_type = std::ptrdiff_t;

    explicit Entry(const particle_type &particle)
            : pos(particle.getPos()), force(force_type()), type(particle.getType()), deactivated(false),
              displacement(0), id(particle.getId()) {}

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

class SCPUParticleData {
public:

    using entry_type = Entry;
    using entries_vec = std::vector<Entry>;
    using entry_index = entries_vec::size_type;
    using new_entries = std::vector<Entry>;
    using particle_type = readdy::model::Particle;
    using top_particle_type = readdy::model::TopologyParticle;
    using force = particle_type::pos_type;
    using displacement = scalar;
    using iterator = entries_vec::iterator;
    using const_iterator = entries_vec::const_iterator;
    using entries_update = std::pair<new_entries, std::vector<entry_index>>;

    SCPUParticleData() = default;
    SCPUParticleData(const SCPUParticleData&) = delete;
    SCPUParticleData& operator=(const SCPUParticleData&) = delete;
    SCPUParticleData(SCPUParticleData&&) = default;
    SCPUParticleData& operator=(SCPUParticleData&&) = default;
    ~SCPUParticleData() = default;

    readdy::model::Particle getParticle(entry_index index) const;

    readdy::model::Particle toParticle(const Entry &e) const {
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
                entries.at(idx) = Entry{p};
            } else {
                entries.emplace_back(p);
            }
        }
    }

    std::vector<entries_vec::size_type> addTopologyParticles(const std::vector<top_particle_type> &particles);

    void removeParticle(const particle_type &particle);

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

    Entry &entry_at(entry_index idx) {
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

    const Entry &entry_at(entry_index idx) const {
        return entries.at(idx);
    }

    const Entry &centry_at(entry_index idx) const {
        return entries.at(idx);
    }

    entry_index addEntry(Entry entry) {
        if(!_blanks.empty()) {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            entries.at(idx) = entry;
            return idx;
        }
        entries.push_back(entry);
        return entries.size()-1;
    }

    void removeEntry(entry_index entry);

    std::vector<entry_index> update(entries_update&&);
    
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
