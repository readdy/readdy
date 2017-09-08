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
 *
 *
 * @file SubCell.cpp
 * @brief 
 * @author clonker
 * @date 4/21/17
 */

#include <readdy/kernel/cpu/nl/SubCell.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {


const bool SubCell::is_leaf() const {
    return _is_leaf;
}

SubCell::SubCell(CellContainer *const super_cell, const vec3 &offset)
        : super(super_cell->data(), super_cell->context(), super_cell->config()) {
    super::_level = static_cast<level_t>(super_cell->level() + 1);
    super::_super_cell = super_cell;
    super::_offset = offset;
}

void SubCell::update_displacements() {
    _maximal_displacements[0] = 0;
    _maximal_displacements[1] = 0;
    if (!_is_leaf) {
        for (auto &cell : _sub_cells) {
            cell.update_displacements();
            const auto &updated_displacements = cell.maximal_displacements();
            if (updated_displacements[0] > _maximal_displacements[0]) {
                _maximal_displacements[0] = updated_displacements[0];
                if(updated_displacements[1] > _maximal_displacements[1]) {
                    _maximal_displacements[1] = updated_displacements[1];
                }
            } else if (updated_displacements[0] > _maximal_displacements[1]) {
                _maximal_displacements[1] = updated_displacements[0];
            }
        }
    } else {
        for (const auto particle_index : _particles_list) {
            const auto &entry = data().entry_at(particle_index);
            assert(!entry.is_deactivated());
            if (entry.displacement > _maximal_displacements[0]) {
                _maximal_displacements[0] = entry.displacement;
            } else if (entry.displacement > _maximal_displacements[1]) {
                _maximal_displacements[1] = entry.displacement;
            }
        }
    }
}

void SubCell::subdivide(const scalar desired_cell_width) {
    throw std::logic_error("this method should not be called on SubCell");
}

void SubCell::refine_uniformly() {
    if (is_leaf()) {
        _is_leaf = false;
        for (std::uint8_t i = 0; i < 3; ++i) {
            _n_sub_cells[i] = 2;
            _sub_cell_size[i] = _size[i] / static_cast<scalar>(_n_sub_cells[i]);
        }
        _sub_cells.reserve(_n_sub_cells[0] * _n_sub_cells[1] * _n_sub_cells[2]);
        for (cell_index i = 0; i < _n_sub_cells[0]; ++i) {
            for (cell_index j = 0; j < _n_sub_cells[1]; ++j) {
                for (cell_index k = 0; k < _n_sub_cells[2]; ++k) {
                    _sub_cells.emplace_back(this,
                                            vec3{_offset.x + i * _sub_cell_size.x, _offset.y + j * _sub_cell_size.y,
                                                 _offset.z + k * _sub_cell_size.z});
                    _sub_cells.back()._size = _sub_cell_size;
                    _sub_cells.back()._contiguous_index = get_contiguous_index(i, j, k, _n_sub_cells[1],
                                                                               _n_sub_cells[2]);
                }
            }
        }
    } else {
        for (auto &sub_cell : _sub_cells) {
            sub_cell.refine_uniformly();
        }
    }
}

void SubCell::update_sub_cell_displacements() {}

void SubCell::setup_uniform_neighbors(const std::uint8_t radius) {
    if (!is_leaf()) {
        std::for_each(_sub_cells.begin(), _sub_cells.end(), [=](sub_cell &sub_cell) {
            sub_cell.setup_uniform_neighbors(static_cast<std::uint8_t>(2 * _level));
        });
    }
    {
        const auto &pbc = _context.applyPBCFun();
        // center w.r.t. simulation coordinates
        auto my_global_center = _offset + .5 * _size - .5 * _root_size;
        auto my_root = root();
        for (int i = -radius; i <= radius; ++i) {
            //auto shifted_pos_i = pbc(my_global_center + vec3{i * _size.x, 0, 0});
            for (int j = -radius; j <= radius; ++j) {
                //auto shifted_pos_j = pbc(shifted_pos_i + vec3{0, j * _size.y, 0});
                for (int k = -radius; k <= radius; ++k) {
                    if (!(i == 0 && j == 0 && k == 0)) {
                        const auto shifted_pos = pbc({my_global_center.x + i * _size.x,
                                                      my_global_center.y + j * _size.y,
                                                      my_global_center.z + k * _size.z});
                        auto cell = my_root->sub_cell_for_position(shifted_pos, _level);
                        if (cell != nullptr && cell != this) {
                            _neighbors.push_back(cell);
                        }
                    }
                }
            }
        }
    }
}

void SubCell::insert_particle(const CellContainer::particle_index index, bool mark_dirty) const {
    _particles_list.add(index);
    if(mark_dirty) set_dirty();
}

void SubCell::insert_particle(const CellContainer::particle_index index, bool mark_dirty) {
    _particles_list.add(index);
    if(mark_dirty) set_dirty();
}

void SubCell::clear() {
    if (!_is_leaf) {
        std::for_each(_sub_cells.begin(), _sub_cells.end(), [](sub_cell &cell) {
            cell.clear();
        });
    } else {
        _particles_list.clear();
    }
}

const ParticlesList &SubCell::particles() const {
    return _particles_list;
}

ParticlesList::particle_indices SubCell::collect_contained_particles() const {
    if (!is_leaf()) {
        ParticlesList::particle_indices result;
        std::vector<ParticlesList::particle_indices> sub_indices;
        sub_indices.reserve(n_sub_cells_total());
        std::size_t total_size = 0;
        for (const auto &sub_cell : _sub_cells) {
            sub_indices.push_back(sub_cell.collect_contained_particles());
            total_size += sub_indices.back().size();
        }
        result.reserve(total_size);
        for (auto &&sub_particles : sub_indices) {
            result.insert(result.end(), std::make_move_iterator(sub_particles.begin()),
                          std::make_move_iterator(sub_particles.end()));
        }
        return result;
    }
    return _particles_list.data();
}

const bool SubCell::neighbor_dirty() const {
    bool neighbor_dirty = false;
    std::for_each(_neighbors.begin(), _neighbors.end(), [&neighbor_dirty](const SubCell *const cell) {
        neighbor_dirty |= cell->is_dirty();
    });
    return neighbor_dirty;
}

const bool SubCell::is_dirty() const {
    return _dirty_flag.get();
}

void SubCell::reset_max_displacements() {
    maximal_displacements()[0] = 0;
    maximal_displacements()[1] = 0;
    for (auto &cell : _sub_cells) {
        cell.reset_max_displacements();
    }
}

void SubCell::reset_particles_displacements() {
    if(!is_leaf()) {
        for(auto& sub_cell : _sub_cells) {
            sub_cell.reset_particles_displacements();
        }
    } else {
        for (const auto p_idx : _particles_list) {
            _data.entry_at(p_idx).displacement = 0;
        }
    }
}

void SubCell::set_dirty() const {
    _dirty_flag.set();
}

void SubCell::unset_dirty() const {
    _dirty_flag.unset();
}

detail::DirtyFlag::DirtyFlag(detail::DirtyFlag &&rhs) noexcept : _is_dirty(rhs._is_dirty.load()) {}

detail::DirtyFlag &detail::DirtyFlag::operator=(detail::DirtyFlag &&rhs) noexcept {
    _is_dirty = rhs.get();
    return *this;
}

bool detail::DirtyFlag::get() const {
    return _is_dirty.load();
}

void detail::DirtyFlag::set() const {
    std::atomic_exchange(&_is_dirty, true);
}

void detail::DirtyFlag::unset() const {
    std::atomic_exchange(&_is_dirty, false);
}
}
}
}
}