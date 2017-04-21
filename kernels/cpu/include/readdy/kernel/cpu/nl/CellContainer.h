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
 * @file CellContainer.h
 * @brief 
 * @author clonker
 * @date 4/21/17
 */
#pragma once

#include <vector>
#include <array>

#include <readdy/model/Vec3.h>
#include <readdy/model/KernelContext.h>
#include <readdy/common/thread/Config.h>

namespace readdy {
namespace kernel {
namespace cpu {

namespace model {
class CPUParticleData;
}


namespace nl {

class SubCell;

class CellContainer {
public:
    using sub_cell = SubCell;
    using sub_cells = std::vector<sub_cell>;
    using cell_index = sub_cells::size_type;
    using dimension = std::array<cell_index, 3>;
    using displacement = scalar;
    using displacement_arr = std::array<displacement, 2>;
    using cell_ref = CellContainer*;
    using cell_ref_list = std::vector<cell_ref>;
    using grid_index = std::tuple<int, int, int>;
    using vec3 = readdy::model::Vec3;

    CellContainer(const model::CPUParticleData& data, const readdy::model::KernelContext& context,
                  const readdy::util::thread::Config& config);

    virtual ~CellContainer() = default;

    CellContainer(CellContainer &&) = default;

    CellContainer &operator=(CellContainer &&) = default;

    CellContainer(const CellContainer &) = delete;

    CellContainer &operator=(const CellContainer &) = delete;

    sub_cells& cells();

    const sub_cells &cells() const;

    dimension &n_cells();

    const dimension &n_cells() const;

    vec3& size();

    const vec3& size() const;

    cell_ref_list& neighbors();

    const cell_ref_list &neighbors() const;

    displacement_arr &maximal_displacements();

    const displacement_arr &maximal_displacements() const;

    cell_index &contiguous_index();

    const cell_index &contiguous_index() const;

    /**
     * return the sub cell that corresponds to the position if it is within bounds (in case of no pbc), otherwise null
     * @todo subtract offset of position so that the index becomes relative
     * @param pos the position
     * @return the corresponding sub cell, otherwise null
     */
    sub_cell* const cell_for_position(const vec3& pos);

    /**
     * return the sub cell that corresponds to the position if it is within bounds (in case of no pbc), otherwise null
     * @param pos the position
     * @return the corresponding sub cell, otherwise null
     */
    const sub_cell* const cell_for_position(const vec3& pos) const;

    /**
     * Returns the cell based on the contiguous index, see cell_for_position().
     * @param index the index (triple)
     * @return the cell if found, otherwise null
     */
    sub_cell* const cell_for_grid_index(const grid_index& index);

    const sub_cell* const cell_for_grid_index(const grid_index& index) const;

    /**
     * Recursively update displacements by calling update_displacements on the sub cells which in turn do this to their
     * sub cells until we reached a LEAF which does the actual computation, anything else is aggregation.
     * Can be parallelized.
     */
    virtual void update_displacements();

    /**
     * Update the maximal displacements of all cells in the container in parallel
     */
    virtual void update_sub_cell_displacements();

    /**
     * sets up subcells with an edge size of >= desired_cell_width
     * @param desired_cell_width
     */
    virtual void subdivide(const scalar desired_cell_width);

    const model::CPUParticleData& data() const;

    const readdy::model::KernelContext& context() const;

    const readdy::util::thread::Config &config() const;

protected:

    template<typename T, typename Dims>
    static T get_contiguous_index(T i, T j, T k, Dims I, Dims J) {
        return k + j * J + i * I * J;
    };

    sub_cells _cells;
    dimension _n_cells{{0, 0, 0}};
    vec3 _cell_size {0, 0, 0};
    vec3 _size{0, 0, 0};

    cell_ref_list _neighbors{};
    displacement_arr _maximal_displacements {{0, 0}};
    cell_index _contiguous_index;

    // the offset of this cell (lower left corner) w.r.t. the simulation box
    vec3 _offset {0, 0, 0};

    const model::CPUParticleData& _data;
    const readdy::model::KernelContext& _context;
    const readdy::util::thread::Config _config;
};

}
}
}
}