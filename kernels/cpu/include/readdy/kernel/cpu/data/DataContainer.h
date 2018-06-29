/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file DataContainer.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <functional>
#include <readdy/model/Context.h>
#include <readdy/common/thread/Config.h>
#include <readdy/common/signals.h>
#include <readdy/common/Utils.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace data {

template<typename T>
class DataContainer {
public:

    using Particle = readdy::model::Particle;
    using TopologyParticle = readdy::model::TopologyParticle;
    using Entries = std::vector<T>;
    using EntriesUpdate = std::vector<T>;
    // tuple of new entries and indices of deleted entries
    using DataUpdate = std::tuple<Entries, std::vector<std::size_t>>;
    using ReorderSignal = readdy::signals::signal<void(const std::vector<std::size_t>)>;
    using topology_index_t = std::ptrdiff_t;
    using size_type = typename Entries::size_type;

    using iterator = typename Entries::iterator;
    using const_iterator = typename Entries::const_iterator;

    DataContainer(const readdy::model::Context &context, thread_pool &pool)
            : _context(context), _pool(pool), reorderSignal(std::make_shared<ReorderSignal>()) {};

    virtual ~DataContainer() = default;

    std::size_t size() const {
        return _entries.size();
    };

    virtual void reserve(std::size_t n) = 0;

    bool empty() const {
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
    std::array<unsigned long long, 3> project(Vec3 value, scalar grid_width) const {
        const auto &box_size = _context.get().boxSize();
        const auto i = static_cast<const unsigned long long>((value.x + .5 * box_size[0]) / grid_width);
        const auto j = static_cast<const unsigned long long>((value.y + .5 * box_size[1]) / grid_width);
        const auto k = static_cast<const unsigned long long>((value.z + .5 * box_size[2]) / grid_width);
        return {{i, j, k}};
    };

    virtual size_type addEntry(T &&entry) = 0;

    virtual void addParticles(const std::vector<Particle> &particles) = 0;

    virtual std::vector<size_type> addTopologyParticles(const std::vector<TopologyParticle> &topologyParticles) = 0;

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
            if(!it_entries->deactivated && it_entries->id == particle.getId()) {
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

    size_type getIndexForId(Particle::id_type id) const {
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

    const Particle::pos_type &pos(size_type index) const {
        return _entries.at(index).pos;
    };

    size_type getNDeactivated() const {
        return _blanks.size();
    }

    readdy::signals::scoped_connection registerReorderEventListener(const ReorderSignal::slot_type &slot) {
        return reorderSignal->connect_scoped(slot);
    }

    virtual std::vector<size_type> update(DataUpdate &&) = 0;

    virtual void displace(size_type entry, const Particle::pos_type &delta) = 0;

    void blanks_moved_to_end() {
        auto n_blanks = _blanks.size();
        std::iota(_blanks.begin(), _blanks.end(), size() - n_blanks - 1);
    }

    void blanks_moved_to_front() {
        std::iota(_blanks.begin(), _blanks.end(), 0);
    }

    const readdy::model::Context &context() const {
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

    std::shared_ptr<ReorderSignal> reorderSignal;
};

struct Entry {
    using Particle = readdy::model::Particle;

    explicit Entry(const Particle &particle)
            : pos(particle.getPos()), force(), type(particle.getType()), deactivated(false), id(particle.getId()) {}

    Entry(Particle::pos_type pos, ParticleTypeId type, Particle::id_type id)
            : pos(pos), type(type), id(id), deactivated(false) {}

    Entry(const Entry &) = default;

    Entry &operator=(const Entry &) = default;

    Entry(Entry &&) noexcept = default;

    Entry &operator=(Entry &&) noexcept = default;

    virtual ~Entry() = default;

    Vec3 force;
    Vec3 pos;
    std::ptrdiff_t topology_index{-1};
    Particle::id_type id;
    Particle::type_type type;
    bool deactivated;
};

using EntryDataContainer = DataContainer<Entry>;

}
}
}
}
