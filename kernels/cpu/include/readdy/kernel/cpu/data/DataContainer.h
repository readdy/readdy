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
#include <readdy/model/KernelContext.h>
#include <readdy/common/thread/Config.h>
#include <readdy/common/signals.h>
#include <readdy/kernel/cpu/util/hilbert.h>
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

    DataContainer(const readdy::model::KernelContext &context, const readdy::util::thread::Config &threadConfig)
            : _context(context), _threadConfig(threadConfig), reorderSignal(std::make_unique<ReorderSignal>()) {};

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

    void hilbertSort(scalar gridWidth) {
        if(!empty()) {
            using indices_it = std::vector<std::size_t>::iterator;
            std::vector<std::size_t> hilbert_indices;
            std::vector<std::size_t> indices(size());

            const auto grainSize = size() / _threadConfig.get().nThreads();
            std::iota(indices.begin(), indices.end(), 0);
            hilbert_indices.resize(size());

            auto worker = [&](std::size_t, const_iterator begin, const_iterator end, indices_it hilbert_begin) {
                {
                    auto it = begin;
                    auto hilbert_it = hilbert_begin;
                    for (; it != end; ++it, ++hilbert_it) {
                        if (!it->deactivated) {
                            const auto integer_coordinates = project(it->pos, gridWidth);
                            *hilbert_it = 1 + static_cast<std::size_t>(hilbert_c2i(3, CHAR_BIT, integer_coordinates.data()));
                        } else {
                            *hilbert_it = 0;
                        }
                    }
                }
            };
            {
                std::vector<std::function<void(std::size_t)>> executables;
                executables.reserve(_threadConfig.get().nThreads());
                const auto& executor = *_threadConfig.get().executor();
                auto data_it = begin();
                auto hilberts_it = hilbert_indices.begin();
                for (std::size_t i = 0; i < _threadConfig.get().nThreads() - 1; ++i) {
                    executables.push_back(executor.pack(worker, data_it, data_it + grainSize, hilberts_it));
                    data_it += grainSize;
                    hilberts_it += grainSize;
                }
                executables.push_back(executor.pack(worker, data_it, end(), hilberts_it));
                executor.execute_and_wait(std::move(executables));
            }
            {
                std::sort(indices.begin(), indices.end(),
                          [&hilbert_indices](std::size_t i, std::size_t j) -> bool {
                              return hilbert_indices[i] < hilbert_indices[j];
                          });
            }

            blanks_moved_to_front();
            std::vector<std::size_t> inverseIndices(indices.size());
            for(std::size_t i = 0; i < indices.size(); ++i) {
                inverseIndices[indices[i]] = i;
            }
            reorderSignal->fire_signal(inverseIndices);
            readdy::util::collections::reorder_destructive(inverseIndices.begin(), inverseIndices.end(), begin());
        }
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

    virtual void displace(T &entry, const Particle::pos_type &delta) = 0;

    void blanks_moved_to_end() {
        auto n_blanks = _blanks.size();
        std::iota(_blanks.begin(), _blanks.end(), size() - n_blanks - 1);
    }

    void blanks_moved_to_front() {
        std::iota(_blanks.begin(), _blanks.end(), 0);
    }

    const readdy::model::KernelContext &context() const {
        return _context.get();
    }

    const readdy::util::thread::Config &threadConfig() const {
        return _threadConfig.get();
    }

protected:
    std::reference_wrapper<const readdy::model::KernelContext> _context;
    std::reference_wrapper<const readdy::util::thread::Config> _threadConfig;

    std::vector<size_type> _blanks {};
    Entries _entries {};

    std::unique_ptr<ReorderSignal> reorderSignal;
};

}
}
}
}