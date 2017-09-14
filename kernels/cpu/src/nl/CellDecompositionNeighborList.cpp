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


#include <readdy/kernel/cpu/nl/CellDecompositionNeighborList.h>

/**
 * << detailed description >>
 *
 * @file CellDecompositionNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 01.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

CellDecompositionNeighborList::CellDecompositionNeighborList(data_type &data,
                                                             const readdy::model::KernelContext &context,
                                                             const readdy::util::thread::Config &config) :
        NeighborList(data, context, config), _cell_container(data, context, config) {}

void CellDecompositionNeighborList::set_up(const util::PerformanceNode &node) {
    auto tx = node.timeit();
    _max_cutoff = _context.get().calculateMaxCutoff();
    _max_cutoff_skin_squared = (_max_cutoff + _skin) * (_max_cutoff + _skin);
    if (_max_cutoff > 0) {
        _cell_container.update_root_size();
        _cell_container.subdivide(_max_cutoff + _skin);
        //_cell_container.refine_uniformly();
        _cell_container.setup_uniform_neighbors();
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

void CellDecompositionNeighborList::update(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if(!_is_set_up) {
        set_up(node.subnode("set_up"));
    }
    _cell_container.clear();
    {
        auto tt = node.subnode("update CLL").timeit();
        fill_container();
    }
    {
        auto tt = node.subnode("fill verlet list").timeit();
        fill_verlet_list();
    }
}

void CellDecompositionNeighborList::fill_container() {
    if (_max_cutoff > 0) {
        const auto &data = _data.get();
        const auto grainSize = data.size() / _config.get().nThreads();
        const auto data_begin = data.cbegin();
        auto worker = [data_begin](std::size_t, data_type::const_iterator begin,
                          data_type::const_iterator end, const CellContainer &container) {
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

void CellDecompositionNeighborList::fill_verlet_list() {
    using namespace std::placeholders;

    if (_max_cutoff > 0) {
        _cell_container.execute_for_each_sub_cell(std::bind(&CellDecompositionNeighborList::fill_cell_verlet_list, this, _1));
    }
}

void CellDecompositionNeighborList::fill_cell_verlet_list(const CellContainer::sub_cell &cell) {
    const auto &d2 = _context.get().distSquaredFun();
    auto &data = _data.get();
    for (const auto particle_index : cell.particles().data()) {
        auto &neighbors = data.neighbors_at(particle_index);
        auto &entry = data.entry_at(particle_index);
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

void CellDecompositionNeighborList::clear(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if (_max_cutoff > 0) {
        _cell_container.clear();
        for (auto &neighbors : _data.get().neighbors()) {
            neighbors.clear();
        }
    }
}

void CellDecompositionNeighborList::updateData(data_type::DataUpdate &&update) {
    _data.get().update(std::forward<data_type::DataUpdate>(update));
}

bool CellDecompositionNeighborList::is_adaptive() const {
    return false;
}

NeighborList::const_iterator CellDecompositionNeighborList::cbegin() const {
    return NeighborListIterator{_data.get().neighbors().begin()};
}

NeighborList::const_iterator CellDecompositionNeighborList::cend() const {
    return NeighborListIterator{_data.get().neighbors().end()};
}


}
}
}
}
