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
 * @file ParticleData.h
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#pragma once
#include <cstddef>
#include <atomic>
#include <vector>
#include <memory>
#include <stack>

#include <readdy/model/Particle.h>
#include <readdy/model/KernelContext.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace model {

class CPUParticleData {
    friend class CPUNeighborList;
public:

    struct Entry;
    using Neighbor = std::size_t;
    using ctx_t = readdy::model::KernelContext;
    using particle_type = readdy::model::Particle;
    using entries_t = std::vector<Entry>;
    using entries_update_t = std::vector<Entry>;
    using top_particle_type = readdy::model::TopologyParticle;
    using neighbors_t = std::vector<Neighbor>;
    using neighbor_list_t = std::vector<neighbors_t>;
    using index_t = entries_t::size_type;
    using update_t = std::pair<entries_update_t, std::vector<index_t>>;
    using force_t = readdy::model::Vec3;
    using displacement_t = double;

    using iterator = decltype(std::declval<entries_t>().begin());
    using const_iterator = decltype(std::declval<entries_t>().cbegin());

    /**
     * Particle data entry with padding such that it fits exactly into 64 bytes.
     */
    struct Entry {
        Entry(const particle_type &particle) : pos(particle.getPos()), force(force_t()), type(particle.getType()),
                                               deactivated(false), displacement(0), id(particle.getId()) {
        }
        Entry(particle_type::pos_type pos, particle_type_type type, particle_type::id_type id)
                : pos(std::move(pos)), type(std::move(type)), id(std::move(id)), deactivated(false),
                  force(), displacement(0) {}

        Entry(const Entry&) = delete;
        Entry& operator=(const Entry&) = delete;
        Entry(Entry&&);
        Entry& operator=(Entry&&);

        bool is_deactivated() const;
        const particle_type::pos_type &position() const;

        force_t force; // 3*8 = 24 bytes
        displacement_t displacement; // 24 + 8 = 32 bytes

    private:
        friend class readdy::kernel::cpu::model::CPUParticleData;
        friend class CPUNeighborList;

        particle_type::pos_type pos; // 32 + 3*8 = 56 bytes
    public:
        particle_type::id_type id; // 56 + 8 = 64
        particle_type::type_type type; // 56 + 4 = 60 bytes
    private:
        bool deactivated; // 60 + 1 = 61 bytes
        char padding[3] {0,0,0}; // 61 + 3 = 64 bytes
    };
    // ctor / dtor
    CPUParticleData(readdy::model::KernelContext *const context);

    ~CPUParticleData();

    // delete move and copy
    CPUParticleData(CPUParticleData &&rhs) = delete;

    CPUParticleData &operator=(CPUParticleData &&rhs) = delete;

    CPUParticleData(const CPUParticleData &rhs) = delete;

    CPUParticleData &operator=(const CPUParticleData &rhs) = delete;

    std::size_t size() const;

    void reserve(std::size_t n);

    bool empty() const;

    void clear();

    void addParticle(const particle_type &particle);

    index_t addEntry(Entry &&entry);

    void addParticles(const std::vector<particle_type> &particles);

    std::vector<entries_t::size_type> addTopologyParticles(const std::vector<top_particle_type> &particles);

    readdy::model::Particle getParticle(const index_t index) const;

    readdy::model::Particle toParticle(const Entry& e) const;

    void removeParticle(const particle_type &particle);

    void removeParticle(const index_t index);

    void removeEntry(index_t entry);

    index_t getIndexForId(const particle_type::id_type) const;

    iterator begin();
    iterator end();
    const_iterator cbegin() const;
    const_iterator cend() const;
    const_iterator begin() const;
    const_iterator end() const;

    Entry& entry_at(index_t);
    const Entry& entry_at(index_t) const;
    const Entry& centry_at(index_t) const;

    neighbors_t& neighbors_at(index_t);
    const neighbors_t& neighbors_at(index_t) const;
    const neighbors_t& cneighbors_at(index_t) const;

    const particle_type::pos_type& pos(index_t) const;

    void setFixPosFun(const ctx_t::fix_pos_fun&);

    void setPBCFun(const ctx_t::pbc_fun& f);

    const ctx_t::fix_pos_fun &fixPosFun() const;

    const ctx_t::pbc_fun &pbcFun() const;

    index_t getNDeactivated() const;

    /**
     *
     * @return vector of new entries
     */
    std::vector<index_t> update(update_t&&);
    void displace(Entry&, const particle_type::pos_type& delta);
    void blanks_moved_to_end();
    void blanks_moved_to_front();

protected:

    std::vector<index_t> blanks;
    neighbor_list_t neighbors;
    entries_t entries;
    ctx_t::fix_pos_fun fixPos;
    ctx_t::pbc_fun pbc;
};

}
}
}
}
