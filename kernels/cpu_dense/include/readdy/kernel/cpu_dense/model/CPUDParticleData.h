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
 * @date 23.11.16
 */

#pragma once
#include <cstddef>
#include <atomic>
#include <vector>
#include <memory>
#include <readdy/model/Particle.h>
#include <readdy/model/KernelContext.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace model {

class CPUDParticleData {
public:

    using particle_type = readdy::model::Particle;

    /**
    * Particle data entry with padding such that it fits exactly into 64 bytes.
    */
    struct Entry {
        Entry(const particle_type &particle) : id(particle.getId()), pos(particle.getPos()),
                                               force(readdy::model::Vec3()), type(particle.getType()),
                                               deactivated(false) { }

        particle_type::id_type id; // 4 bytes
        particle_type::pos_type pos; // 4 + 24 = 28 bytes
        readdy::model::Vec3 force; // 28 + 24 = 52 bytes
        particle_type::type_type type; // 52 + 8 = 60 bytes

        const particle_type::pos_type& position() const {
            return pos;
        }

        bool deactivated; // 60 + 1 = 61 bytes
    private:
        friend class readdy::kernel::cpu_dense::model::CPUDParticleData;

        bool padding[3]; // 61 + 3 = 64 bytes
    };

    using marked_count_t = std::atomic<std::size_t>;
    using entries_t = std::vector<Entry>;
    using index_t = entries_t::size_type;
    using iterator = entries_t::iterator;
    using const_iterator = entries_t::const_iterator;
    using update_t = std::vector<Entry>;

    // ctor / dtor
    CPUDParticleData(const readdy::model::KernelContext*const);

    CPUDParticleData(const readdy::model::KernelContext*const, unsigned int capacity);

    ~CPUDParticleData();

    // move
    CPUDParticleData(CPUDParticleData &&rhs) = delete;

    CPUDParticleData &operator=(CPUDParticleData &&rhs) = delete;

    // copy
    CPUDParticleData(const CPUDParticleData &rhs) = delete;

    CPUDParticleData &operator=(const CPUDParticleData &rhs) = delete;

    iterator begin();

    const_iterator begin() const;

    const_iterator cbegin() const;

    iterator end();

    const_iterator end() const;

    const_iterator cend() const;

    std::size_t size() const;

    std::size_t max_size() const;

    Entry& entry_at(const index_t);
    const Entry& entry_at(const index_t) const;
    const Entry& centry_at(const index_t) const;

    readdy::model::Particle toParticle(const Entry& e) const;

    bool empty() const;

    void clear();

    void addParticle(const particle_type &particle);

    void addParticles(const std::vector<particle_type> &particles);

    void displace(Entry&, const particle_type::pos_type& delta);

    void update(update_t&&);

    /**
     * Remove a particle via its unique id.
     * @param particle the particle to be removed
     */
    void removeParticle(const particle_type &particle);

    void removeParticle(const std::size_t index);

    bool isMarkedForDeactivation(const std::size_t index);

    std::size_t getDeactivatedIndex() const;

    std::size_t getNDeactivated() const;

    /**
     * This method is the counterpart to markForDeactivation.
     * The particles that were marked are now deactivated, i.e.,
     * for each marked particle:
     *   - If it is at the very end of the particle list, the
     *     counters are updated.
     *   - If not, the particle is swapped with the last active particle,
     *     so that again, all deactivated particles reside at the end
     *     of the internal data structure.
     */
    void deactivateMarked();

    void deactivate(index_t);

    void deactivate(Entry&);

protected:
    entries_t entries;
    size_t deactivated_index;
    size_t n_deactivated;
    const readdy::model::KernelContext*const ctx;
    mutable marked_count_t n_marked;
};


}
}
}
}
