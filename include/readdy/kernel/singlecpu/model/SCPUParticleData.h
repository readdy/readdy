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

    Entry(const Entry &) = delete;

    Entry &operator=(const Entry &) = delete;

    Entry(Entry &&) noexcept = default;

    Entry &operator=(Entry &&) noexcept = default;

    ~Entry() = default;

    bool is_deactivated() const;

    const particle_type::pos_type &position() const;

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

    readdy::model::Particle toParticle(const Entry &e) const;

    void addParticle(const particle_type &particle);

    void addParticles(const std::vector<particle_type> &particles);

    std::vector<entries_vec::size_type> addTopologyParticles(const std::vector<top_particle_type> &particles);

    void removeParticle(const particle_type &particle);

    void removeParticle(size_t index);

    iterator begin();

    iterator end();

    const_iterator cbegin() const;

    const_iterator cend() const;

    const_iterator begin() const;

    const_iterator end() const;

    Entry &entry_at(entry_index);

    entry_index size() const;

    entry_index n_deactivated() const;

    void reserve(std::size_t n);

    void clear();

    const Entry &entry_at(entry_index) const;

    const Entry &centry_at(entry_index) const;

    entry_index addEntry(Entry &&entry);

    void removeEntry(entry_index entry);

    std::vector<entry_index> update(entries_update&&);

protected:

    entries_vec entries;
    std::vector<entry_index> blanks;

};


}
}
}
}
