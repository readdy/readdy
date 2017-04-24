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
    if (!_is_leaf) {
        _maximal_displacements[0] = 0;
        _maximal_displacements[1] = 0;
        for (auto &cell : _sub_cells) {
            cell.update_displacements();
            const auto &updated_displacements = cell.maximal_displacements();
            if (updated_displacements[0] > _maximal_displacements[0]) {
                _maximal_displacements[0] = updated_displacements[0];
                _maximal_displacements[1] = std::max(_maximal_displacements[1], updated_displacements[1]);
            } else if (updated_displacements[0] > _maximal_displacements[1]) {
                _maximal_displacements[1] = updated_displacements[0];
            }
        }
    } else {
        // todo update displacements wrt. contained particles
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
    } else {
        const auto &pbc = _context.getPBCFun();
        // center w.r.t. simulation coordinates
        auto my_global_center = _offset + .5 * _size - .5 * _root_size;
        auto my_root = root();
        for (int i = -radius; i <= radius; ++i) {
            auto shifted_pos_i = pbc(my_global_center + vec3{i * _size.x, 0, 0});
            for (int j = -radius; j <= radius; ++j) {
                auto shifted_pos_j = pbc(shifted_pos_i + vec3{0, j * _size.y, 0});
                for (int k = -radius; k <= radius; ++k) {
                    if (!(i == 0 && j == 0 && k == 0)) {
                        //const auto shifted_pos = pbc(shifted_pos_j + vec3{0, 0, k*_size.z});
                        //if(shifted_pos == shifted_pos_j) break;
                        const auto shifted_pos = pbc({my_global_center.x + i * _size.x,
                                                      my_global_center.y + j * _size.y,
                                                      my_global_center.z + k * _size.z});
                        auto cell = my_root->leaf_cell_for_position(shifted_pos);
                        if (cell && cell != this) {
                            _neighbors.push_back(cell);
                        }
                    }
                }
            }
        }
    }
}

}
}
}
}