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
 * @file NeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 24.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/cpu_legacy/nl/AdaptiveNeighborList.h>

namespace readdy {
namespace kernel {
namespace cpu_legacy {
namespace nl {

AdaptiveNeighborList::AdaptiveNeighborList(data::EntryDataContainer *data,
                                           const readdy::model::Context &context,
                                           const readdy::util::thread::Config &config, bool hilbert_sort)
        : NeighborList(context, config), _data(data), _cell_container(_data, context, config),
          _hilbert_sort(hilbert_sort) {}

AdaptiveNeighborList::AdaptiveNeighborList(const readdy::model::Context &context,
                           const readdy::util::thread::Config &config, bool hilbert_sort)
        : NeighborList(context, config), _data(context, config), _cell_container(_data, context, config),
          _hilbert_sort(hilbert_sort) {}

const CellContainer &AdaptiveNeighborList::cell_container() const {
    return _cell_container;
}

void AdaptiveNeighborList::set_up(const util::PerformanceNode &node) {
    auto tx = node.timeit();
    _max_cutoff = _context.get().calculateMaxCutoff();
    _max_cutoff_skin_squared = (_max_cutoff + _skin) * (_max_cutoff + _skin);
    if (_max_cutoff > 0) {
        _cell_container.update_root_size();
        _cell_container.subdivide(_max_cutoff + _skin);
        //_cell_container.refine_uniformly();
        _cell_container.setup_uniform_neighbors();
        if (_hilbert_sort) {
            _data.hilbertSort(_max_cutoff + _skin);
        }
        {
            auto t = node.subnode("fill_container").timeit();
            fill_container();
        }
        {
            auto t = node.subnode("fill_verlet_list").timeit();
            fill_verlet_list();
        }   
    }
    _is_set_up = true;
}

void AdaptiveNeighborList::update(const util::PerformanceNode &node) {
    if(!_is_set_up) {
        set_up(node.subnode("set_up"));
    } else {
        if (_max_cutoff > 0) {
            bool too_far = !_cell_container.update_sub_cell_displacements_and_mark_dirty(_max_cutoff, _skin);
            bool too_many = _cell_container.n_dirty_macro_cells() >= .9 * _cell_container.n_sub_cells_total();
            log::trace("updating {}% ({} of {})",
                       100. * _cell_container.n_dirty_macro_cells() / _cell_container.n_sub_cells_total(),
                       _cell_container.n_dirty_macro_cells(), _cell_container.n_sub_cells_total());
            if (!too_far && !too_many) {
                _cell_container.update_dirty_cells();
                handle_dirty_cells();
            } else {
                if (too_far) {
                    // a particle travelled too far, we need to re-setup the whole thing
                    log::warn("A particle's displacement has been more than r_c + r_s = {} + {} = {}, which means that "
                                      "it might have left its cell linked-list cell. This should, if at all, only happen "
                                      "very rarely and triggers a complete rebuild of the neighbor list.",
                              _max_cutoff, _skin, _max_cutoff + _skin);
                }
                if (too_many) {
                    log::debug(
                            "More than 90% of the cells were marked dirty, thus re-create the whole neighbor list rather"
                                    " than update it adaptively");
                }
                clear_cells();
                fill_container();
                fill_verlet_list();
            }
        }
    }
}

const readdy::util::thread::Config &AdaptiveNeighborList::config() const {
    return _config;
}

void AdaptiveNeighborList::fill_container() {
    if (_max_cutoff > 0) {
        const auto &data = _data;
        const auto grainSize = data.size() / _config.get().nThreads();
        const auto data_begin = data.cbegin();
        auto worker = [data_begin](std::size_t, data::EntryDataContainer::const_iterator begin,
                                   data::EntryDataContainer::const_iterator end, const CellContainer &container) {
            auto i = std::distance(data_begin, begin);
            for (auto it = begin; it != end; ++it, ++i) {
                if (!it->deactivated) {
                    container.insert_particle(static_cast<const CellContainer::particle_index>(i));
                }
            }
        };

        auto& executor = *_config.get().executor();

        std::vector<std::function<void(std::size_t)>> executables;
        executables.reserve(_config.get().nThreads());

        auto it = data.cbegin();
        for (auto i = 0; i < _config.get().nThreads() - 1; ++i) {
            executables.push_back(executor.pack(worker, it, it + grainSize, std::cref(_cell_container)));
            it += grainSize;
        }
        executables.push_back(executor.pack(worker, it, data.end(), std::cref(_cell_container)));

        executor.execute_and_wait(std::move(executables));
    }
}

void AdaptiveNeighborList::clear_cells() {
    _cell_container.clear();
}

void AdaptiveNeighborList::clear(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if (_max_cutoff > 0) {
        clear_cells();
        for (auto &neighbors : _data.neighbors()) {
            neighbors.clear();
        }
    }
}

void AdaptiveNeighborList::fill_verlet_list() {
    using namespace std::placeholders;
    if (_max_cutoff > 0) {
        _cell_container.execute_for_each_sub_cell(std::bind(&AdaptiveNeighborList::fill_cell_verlet_list, this, _1, true));
    }
}

void AdaptiveNeighborList::fill_cell_verlet_list(const CellContainer::sub_cell &cell, const bool reset_displacement) {

    const auto &d2 = _context.get().distSquaredFun();
    auto &data = _data;
    for (const auto particle_index : cell.particles().data()) {
        auto &neighbors = data.neighbors_at(particle_index);
        auto &entry = data.entry_at(particle_index);
        if (reset_displacement) {
            data.displacements().at(particle_index) = 0;
        }
        neighbors.clear();
        for (const auto p_i : cell.particles().data()) {
            if (p_i != particle_index) {
                const auto distSquared = d2(entry.pos, data.pos(p_i));
                if (distSquared < _max_cutoff_skin_squared) {
                    neighbors.push_back(p_i);
                }
            }
        }
        for (const auto &neighbor_cell : cell.neighbors()) {
            for (const auto p_j : neighbor_cell->particles().data()) {
                const auto distSquared = d2(entry.pos, data.pos(p_j));
                if (distSquared < _max_cutoff_skin_squared) {
                    neighbors.push_back(p_j);
                }
            }
        }
    }
}

bool &AdaptiveNeighborList::performs_hilbert_sort() {
    return _hilbert_sort;
}

const bool &AdaptiveNeighborList::performs_hilbert_sort() const {
    return _hilbert_sort;
}

void AdaptiveNeighborList::sort_by_hilbert_curve() {
    clear({});
    if (_hilbert_sort) {
        _data.hilbertSort(_max_cutoff + _skin);
    } else {
        log::error("requested hilbert sort but it is turned off!");
    }
    fill_container();
    fill_verlet_list();
}

void AdaptiveNeighborList::updateData(data_type::DataUpdate &&update) {
    const auto& decayed_particles = std::get<1>(update);
    for(const auto p_idx : decayed_particles) {
        const auto sub_cell = _cell_container.leaf_cell_for_position(_data.pos(p_idx));
        if (sub_cell != nullptr) {
            const auto& particles = sub_cell->particles();
            if(particles.erase_if_found(p_idx)) {
                sub_cell->super_cell()->set_dirty();
            } else {
                bool found = false;
                for(const auto neighbor : sub_cell->neighbors()) {
                    if(!found) {
                        const auto& neighbor_particles = neighbor->particles();
                        if(neighbor_particles.erase_if_found(p_idx)) {
                            neighbor->super_cell()->set_dirty();
                            found = true;
                        }
                    } else {
                        break;
                    }
                }
                if(!found) {
                    log::critical("something went very very wrong");
                }
            }
        }
    }

    auto new_entries = _data.update(std::move(update));

    for(const auto p_idx : new_entries) {
        _cell_container.insert_particle(p_idx, true);
    }

    handle_dirty_cells();
}

void AdaptiveNeighborList::handle_dirty_cells() {
    _cell_container.execute_for_each_sub_cell([this](const CellContainer::sub_cell &cell) {
        if (cell.is_dirty() || cell.neighbor_dirty() ) { // or neighbor? ||
            fill_cell_verlet_list(cell, false);
        }
    });
    _cell_container.unset_dirty();
}

bool AdaptiveNeighborList::is_adaptive() const {
    return true;
}

NeighborList::const_iterator AdaptiveNeighborList::cbegin() const {
    return NeighborListIterator{_data.neighbors().begin(), _data.neighbors().end(), true};
}

NeighborList::const_iterator AdaptiveNeighborList::cend() const {
    return NeighborListIterator{_data.neighbors().end(), _data.neighbors().end(), true};
}

data::EntryDataContainer *AdaptiveNeighborList::data() {
    return &_data;
}

const data::EntryDataContainer *AdaptiveNeighborList::data() const {
    return &_data;
}

const data::NLDataContainer &AdaptiveNeighborList::nlData() const {
    return _data;
}

size_t AdaptiveNeighborList::size() const {
    return _data.neighbors().size();
}

}
}
}
}
