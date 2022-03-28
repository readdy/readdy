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
 * @file DataContainer.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.09.17
 * @copyright BSD-3
 */

#pragma once

#include <functional>
#include <readdy/model/Context.h>
#include <readdy/common/thread/Config.h>
#include <readdy/common/signals.h>
#include <readdy/common/Utils.h>

namespace readdy::kernel::cpu::data {

template<typename T>
class DataContainer {
public:

    using Particle = readdy::model::Particle;
    using Entries = std::vector<T>;
    using EntriesUpdate = std::vector<T>;
    // tuple of new entries and indices of deleted entries
    using DataUpdate = std::tuple<Entries, std::vector<std::size_t>>;
    using topology_index_t = std::ptrdiff_t;
    using size_type = typename Entries::size_type;

    using iterator = typename Entries::iterator;
    using const_iterator = typename Entries::const_iterator;

    DataContainer(const readdy::model::Context &context, thread_pool &pool)
            : _context(context), _pool(pool) {};

    virtual ~DataContainer() = default;

    [[nodiscard]] std::size_t size() const {
        return _entries.size();
    };

    virtual void reserve(std::size_t n) = 0;

    [[nodiscard]] bool empty() const {
        return size() == getNDeactivated();
    };

    void clear() {
        _entries.clear();
        _blanks.clear();
    };

    void addParticle(const Particle &particle) {
        addParticles({particle});
    };

    /**
     * project into an unsigned long long int assuming that value is within the sim box
     * @param value the value
     * @param grid_width precision
     * @return bitmask
     */
    [[nodiscard]] std::array<unsigned long long, 3> project(Vec3 value, scalar grid_width) const {
        const auto &box_size = _context.get().boxSize();
        const auto i = static_cast<const unsigned long long>((value.x + .5 * box_size[0]) / grid_width);
        const auto j = static_cast<const unsigned long long>((value.y + .5 * box_size[1]) / grid_width);
        const auto k = static_cast<const unsigned long long>((value.z + .5 * box_size[2]) / grid_width);
        return {{i, j, k}};
    };

    virtual size_type addEntry(T &&entry) = 0;

    virtual void addParticles(const std::vector<Particle> &particles) = 0;

    virtual std::vector<size_type> addTopologyParticles(const std::vector<Particle> &topologyParticles) = 0;

    Particle getParticle(size_type index) const {
        const auto& entry = *(_entries.begin() + index);
        if(entry.deactivated) {
            log::error("Requested deactivated particle at index {}!", index);
        }
        return toParticle(entry);
    };

    Particle toParticle(const T &e) const {
        return readdy::model::Particle(e.pos, e.type, e.id);
    };

    void removeParticle(const Particle &particle) {
        auto it_entries = begin();
        std::size_t idx = 0;
        for(; it_entries != end(); ++it_entries, ++idx) {
            if(!it_entries->deactivated && it_entries->id == particle.id()) {
                _blanks.push_back(idx);
                it_entries->deactivated = true;
                return;
            }
        }
        log::error("Tried to remove particle ({}) which did not exist or was already deactivated!", particle);
    };

    void removeParticle(size_type index) {
        auto& p = *(_entries.begin() + index);
        if(!p.deactivated) {
            _blanks.push_back(index);
            p.deactivated = true;
        } else {
            log::error("Tried to remove particle (index={}), that was already removed!", index);
        }
    };

    void removeEntry(size_type index) {
        auto &entry = _entries.at(index);
        if(!entry.deactivated) {
            entry.deactivated = true;
            _blanks.push_back(index);
        } else {
            log::critical("Tried removing particle {} which was already deactivated!", index);
        }
    };

    size_type getIndexForId(ParticleId id) const {
        auto find_it = std::find_if(_entries.begin(), _entries.end(), [id](const T& e) {
            return !e.deactivated && e.id == id;
        });
        if(find_it != _entries.end()) {
            return static_cast<size_type>(std::distance(_entries.begin(), find_it));
        }
        throw std::out_of_range("requested id was not to be found in particle data");
    };

    virtual iterator begin() {
        return _entries.begin();
    };

    virtual iterator end() {
        return _entries.end();
    };

    virtual const_iterator cbegin() const {
        return _entries.cbegin();
    };

    virtual const_iterator cend() const {
        return _entries.cend();
    };

    virtual const_iterator begin() const {
        return _entries.cbegin();
    };

    virtual const_iterator end() const {
        return _entries.cend();
    };

    T &entry_at(size_type index) {
        return _entries.at(index);
    }

    const T &entry_at(size_type index) const {
        return _entries.at(index);
    }

    const T &centry_at(size_type index) const {
        return _entries.at(index);
    }

    const Particle::Position &pos(size_type index) const {
        return _entries.at(index).pos;
    };

    size_type getNDeactivated() const {
        return _blanks.size();
    }

    virtual std::vector<size_type> update(DataUpdate &&) = 0;

    virtual void displace(size_type entry, const Particle::Position &delta) = 0;

    [[nodiscard]] const readdy::model::Context &context() const {
        return _context.get();
    }

    auto &pool() {
        return _pool.get();
    }

    const auto &entries() const {
        return _entries;
    }

    const std::vector<size_type> &blanks() const {
        return _blanks;
    }

protected:
    std::reference_wrapper<const readdy::model::Context> _context;
    std::reference_wrapper<thread_pool> _pool;

    std::vector<size_type> _blanks {};
    Entries _entries {};
};

struct Entry {
    using Particle = readdy::model::Particle;

    explicit Entry(const Particle &particle)
            : pos(particle.pos()), force(), type(particle.type()), deactivated(false), id(particle.id()) {}

    Entry(Particle::Position pos, ParticleTypeId type, ParticleId id)
            : pos(pos), type(type), id(id), deactivated(false) {}

    Entry(const Entry &) = default;

    Entry &operator=(const Entry &) = default;

    Entry(Entry &&) noexcept = default;

    Entry &operator=(Entry &&) noexcept = default;

    virtual ~Entry() = default;

    Vec3 force;
    Vec3 pos;
    std::ptrdiff_t topology_index{-1};
    ParticleId id;
    ParticleTypeId type;
    bool deactivated;
};

using EntryDataContainer = DataContainer<Entry>;

}
