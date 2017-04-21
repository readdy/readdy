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
 * @file SubCell.h
 * @brief 
 * @author clonker
 * @date 4/21/17
 */
#pragma once

#include "CellContainer.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class SubCell : public CellContainer {
    using super = CellContainer;
public:
    using particle_ref = int;

    SubCell(const CellContainer &super_cell, const vec3 &offset);

    const bool is_leaf() const;

    virtual void update_displacements() override;

    virtual void subdivide(const scalar desired_cell_width) override;

    void refine_uniformly();

private:
    const CellContainer &super_cell;

    bool _is_leaf {true};

    // change visibility to private, this should not be used with sub cells
    void update_sub_cell_displacements() override;

};

}
}
}
}