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

#include <readdy/common/signals.h>
#include <readdy/common/thread/Config.h>
#include <readdy/model/Particle.h>
#include <readdy/model/KernelContext.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {
class NeighborList;
class AdaptiveNeighborList;
class ContiguousCLLNeighborList;
class DynamicCLLNeighborList;
class CellDecompositionNeighborList;
}
namespace model {

class CPUParticleData {
    friend class CPUNeighborList;
    friend class readdy::kernel::cpu::nl::NeighborList;
    friend class readdy::kernel::cpu::nl::AdaptiveNeighborList;
    friend class readdy::kernel::cpu::nl::ContiguousCLLNeighborList;
    friend class readdy::kernel::cpu::nl::DynamicCLLNeighborList;
    friend class readdy::kernel::cpu::nl::CellDecompositionNeighborList;
public:

    struct Entry;
    using Neighbor = std::size_t;
    using ctx_t = readdy::model::KernelContext;
    using particle_type = readdy::model::Particle;
    using entries_t = std::vector<Entry>;
    using entries_update_t = std::vector<Entry>;
    using entry_t = entries_t::value_type;
    using top_particle_type = readdy::model::TopologyParticle;
    using neighbors_t = std::vector<Neighbor>;
    using neighbor_list_t = std::vector<neighbors_t>;
    using index_t = entries_t::size_type;
    using update_t = std::pair<entries_update_t, std::vector<index_t>>;
    using vec3 = Vec3;
    using force_t = vec3;
    using displacement_t = readdy::scalar;
    using reorder_signal_t = readdy::signals::signal<void(const std::vector<std::size_t>)>;
    using topology_index_t = std::ptrdiff_t;

    using iterator = decltype(std::declval<entries_t>().begin());
    using const_iterator = decltype(std::declval<entries_t>().cbegin());

    using neighbors_list_iterator = neighbor_list_t::iterator;
    using neighbors_list_const_iterator = neighbor_list_t::const_iterator;

    /**
     * Particle data entry with padding such that it fits exactly into 64 bytes.
     */
    struct Entry {
        explicit Entry(const particle_type &particle) : pos(particle.getPos()), force(force_t()), type(particle.getType()),
                                               deactivated(false), displacement(0), id(particle.getId()) {
        }
        Entry(particle_type::pos_type pos, particle_type_type type, particle_type::id_type id)
                : pos(pos), type(type), id(id), deactivated(false), displacement(0) {}

        Entry(const Entry&) = delete;
        Entry& operator=(const Entry&) = delete;
        Entry(Entry&&) noexcept = default;
        Entry& operator=(Entry&&) noexcept = default;
        ~Entry() = default;

        bool is_deactivated() const;
        const particle_type::pos_type &position() const;

        force_t force; // 3*8 = 24 bytes
        displacement_t displacement; // 24 + 8 = 32 bytes

        particle_type::pos_type pos; // 32 + 3*8 = 56 bytes
        topology_index_t topology_index {-1};
        particle_type::id_type id; // 56 + 8 = 64
        particle_type::type_type type; // 56 + 4 = 60 bytes
        bool deactivated; // 60 + 1 = 61 bytes
        //char padding[3] {0,0,0}; // 61 + 3 = 64 bytes
    };
    // ctor / dtor
    CPUParticleData(readdy::model::KernelContext *context, const readdy::util::thread::Config &_config);

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

    /**
     * project into an unsigned long long int assuming that value is within the sim box
     * @param value the value
     * @param grid_width precision
     * @return bitmask
     */
    std::array<unsigned long long, 3> project(vec3 value, scalar grid_width) const;

    void hilbert_sort(scalar grid_width);

    index_t addEntry(Entry &&entry);

    void addParticles(const std::vector<particle_type> &particles);

    std::vector<entries_t::size_type> addTopologyParticles(const std::vector<top_particle_type> &particles);

    readdy::model::Particle getParticle(index_t index) const;

    readdy::model::Particle toParticle(const Entry& e) const;

    void removeParticle(const particle_type &particle);

    void removeParticle(index_t index);

    void removeEntry(index_t entry);

    index_t getIndexForId(particle_type::id_type) const;

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

    index_t getNDeactivated() const;

    readdy::signals::scoped_connection registerReorderEventListener(const reorder_signal_t::slot_type &slot);

    /**
     *
     * @return vector of new entries
     */
    std::vector<index_t> update(update_t&&);
    void displace(Entry&, const particle_type::pos_type& delta);
    void blanks_moved_to_end();
    void blanks_moved_to_front();

    bool& trackDisplacement();

    const bool& trackDisplacement() const;

protected:

    std::vector<index_t> blanks;
    neighbor_list_t neighbors;
    entries_t entries;

    bool _trackDisplacement {true};

    std::unique_ptr<reorder_signal_t> reorderSignal;

    const readdy::model::KernelContext* const context;
    const readdy::util::thread::Config& _config;
};

}
}
}
}
