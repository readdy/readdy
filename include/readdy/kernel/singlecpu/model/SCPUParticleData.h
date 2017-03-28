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
    using entries_t = std::vector<Entry>;
    using particle_type = readdy::model::Particle;
    using force_t = particle_type::pos_type;
    using displacement_t = force_t::value_t;

    Entry(const particle_type &particle) : pos(particle.getPos()), force(force_t()), type(particle.getType()),
                                           deactivated(false), displacement(0), id(particle.getId()) {}

    Entry(const Entry &) = delete;

    Entry &operator=(const Entry &) = delete;

    Entry(Entry &&) = default;

    Entry &operator=(Entry &&) = default;

    bool is_deactivated() const;

    const particle_type::pos_type &position() const;

    force_t force;
    displacement_t displacement;
    particle_type::pos_type pos;

private:
    friend class readdy::kernel::scpu::model::SCPUParticleData;

public:
    particle_type::id_type id;
    particle_type::type_type type;
private:
    bool deactivated;
};

class SCPUParticleData {
public:

    using entries_t = std::vector<Entry>;
    using index_t = entries_t::size_type;
    using entries_update_t = std::vector<Entry>;
    using particle_type = readdy::model::Particle;
    using top_particle_type = readdy::model::TopologyParticle;
    using force_t = particle_type::pos_type;
    using displacement_t = force_t::value_t;
    using iterator = entries_t::iterator;
    using const_iterator = entries_t::const_iterator;
    using update_t = std::pair<entries_update_t, std::vector<index_t>>;

    SCPUParticleData() = default;
    SCPUParticleData(const SCPUParticleData&) = delete;
    SCPUParticleData& operator=(const SCPUParticleData&) = delete;
    SCPUParticleData(SCPUParticleData&&) = default;
    SCPUParticleData& operator=(SCPUParticleData&&) = default;

    readdy::model::Particle getParticle(const index_t index) const;

    readdy::model::Particle toParticle(const Entry &e) const;

    void addParticle(const particle_type &particle);

    void addParticles(const std::vector<particle_type> &particles);

    std::vector<entries_t::size_type> addTopologyParticles(const std::vector<top_particle_type> &particles);

    void removeParticle(const particle_type &particle);

    void removeParticle(const size_t index);

    iterator begin();

    iterator end();

    const_iterator cbegin() const;

    const_iterator cend() const;

    const_iterator begin() const;

    const_iterator end() const;

    Entry &entry_at(index_t);

    index_t size() const;

    index_t n_deactivated() const;

    void reserve(std::size_t n);

    void clear();

    const Entry &entry_at(index_t) const;

    const Entry &centry_at(index_t) const;

    index_t addEntry(Entry &&entry);

    void removeEntry(index_t entry);

    std::vector<index_t> update(update_t&&);

protected:

    entries_t entries;
    std::vector<index_t> blanks;

};


}
}
}
}
