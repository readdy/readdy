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

SubCell::SubCell(const CellContainer &super_cell, const vec3 &offset)
        : super(super_cell.data(), super_cell.context(), super_cell.config()), super_cell(super_cell) {
    super::_offset = offset;
}

void SubCell::update_displacements() {
    if(!_is_leaf) {
        _maximal_displacements[0] = 0;
        _maximal_displacements[1] = 0;
        for(auto& cell : _cells) {
            cell.update_displacements();
            const auto& updated_displacements = cell.maximal_displacements();
            if(updated_displacements[0] > _maximal_displacements[0]) {
                _maximal_displacements[0] = updated_displacements[0];
                _maximal_displacements[1] = std::max(_maximal_displacements[1], updated_displacements[1]);
            } else if(updated_displacements[0] > _maximal_displacements[1]) {
                _maximal_displacements[1] = updated_displacements[0];
            }
        }
    } else {
        // todo update displacements wrt. contained particles
    }
}

void SubCell::subdivide(const scalar desired_cell_width) {
    // todo (set is_leaf to false)
}

void SubCell::refine_uniformly() {
    // todo (set is_leaf to false)
}

void SubCell::update_sub_cell_displacements() {}

}
}
}
}