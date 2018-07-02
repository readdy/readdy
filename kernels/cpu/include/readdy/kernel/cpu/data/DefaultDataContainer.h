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
 * @file DefaultDataContainer.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/model/Particle.h>
#include <readdy/kernel/cpu/pool.h>
#include <readdy/common/boundary_condition_operations.h>
#include "DataContainer.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace data {

class NLDataContainer;

class DefaultDataContainer : public EntryDataContainer {
    using super = DataContainer<readdy::kernel::cpu::data::Entry>;
    friend class NLDataContainer;
public:

    using entry_type = typename super::Entries::value_type;

    explicit DefaultDataContainer(EntryDataContainer *entryDataContainer)
            : DataContainer(entryDataContainer->context(), entryDataContainer->pool()) {
        _entries = entryDataContainer->entries();
        _blanks = entryDataContainer->blanks();
    }

    DefaultDataContainer(const model::Context &context, thread_pool &pool) : DataContainer(context, pool) {};

    void reserve(std::size_t n) override {
        _entries.reserve(n);
    };

    size_type addEntry(Entry &&entry) override {
        if(!_blanks.empty()) {
            const auto idx = _blanks.back();
            _blanks.pop_back();
            _entries.at(idx) = std::move(entry);
            return idx;
        }

        _entries.push_back(std::move(entry));
        return _entries.size()-1;
    }

    void addParticles(const std::vector<Particle> &particles) override {
        for(const auto& p : particles) {
            if(!_blanks.empty()) {
                const auto idx = _blanks.back();
                _blanks.pop_back();
                _entries.at(idx) = Entry(p);
            } else {
                _entries.emplace_back(p);
            }
        }
    }

    std::vector<size_type>
    addTopologyParticles(const std::vector<TopologyParticle> &topologyParticles) override {
        std::vector<size_type> indices;
        indices.reserve(topologyParticles.size());
        for(const auto& p : topologyParticles) {
            if(!_blanks.empty()) {
                const auto idx = _blanks.back();
                _blanks.pop_back();
                _entries.at(idx) = Entry(p);
                indices.push_back(idx);
            } else {
                _entries.emplace_back(p);
                indices.push_back(_entries.size()-1);
            }
        }
        return indices;
    }

    std::vector<size_type> update(DataUpdate &&update) override {
        auto &&newEntries = std::move(std::get<0>(update));
        auto &&removedEntries = std::move(std::get<1>(update));

        auto it_del = removedEntries.begin();
        for(auto&& newEntry : newEntries) {
            if(it_del != removedEntries.end()) {
                _entries.at(*it_del) = std::move(newEntry);
                ++it_del;
            } else {
                addEntry(std::move(newEntry));
            }
        }
        while(it_del != removedEntries.end()) {
            removeEntry(*it_del);
            ++it_del;
        }
        return {};
    }

    void displace(size_type index, const Particle::pos_type &delta) override {
        auto &entry = _entries.at(index);
        entry.pos += delta;
        bcs::fixPosition(entry.pos, _context.get().boxSize().data(),
                         _context.get().periodicBoundaryConditions().data());
    };

    /*void hilbertSort(scalar gridWidth) {
        if(!empty()) {
            using indices_it = std::vector<std::size_t>::iterator;
            std::vector<std::size_t> hilbert_indices;
            std::vector<std::size_t> indices(size());

            const auto grainSize = size() / _pool.get().size();
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
                std::vector<util::thread::joining_future<void>> futures;
                futures.reserve(_pool.get().size());
                auto data_it = begin();
                auto hilberts_it = hilbert_indices.begin();
                for (std::size_t i = 0; i < _pool.get().size() - 1; ++i) {
                    futures.emplace_back(_pool.get().push(worker, data_it, data_it + grainSize, hilberts_it));
                    data_it += grainSize;
                    hilberts_it += grainSize;
                }
                futures.emplace_back(_pool.get().push(worker, data_it, end(), hilberts_it));
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
    };*/

};

}
}
}
}
