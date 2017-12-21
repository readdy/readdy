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

#include <readdy/common/thread/Config.h>
#include <readdy/common/common.h>
#include <readdy/common/Index.h>
#include <readdy/model/Context.h>

#include <readdy/kernel/cpu_legacy/data/NLDataContainer.h>

#include "ParticlesList.h"

namespace readdy {
namespace kernel {
namespace cpu_legacy {
namespace nl {

class SubCell;

class CellContainer {
public:
    using sub_cell = SubCell;
    using sub_cells_t = std::vector<sub_cell>;
    using cell_index = sub_cells_t::size_type;
    using dimension = std::array<cell_index, 3>;
    using displacement = scalar;
    using displacement_arr = std::array<displacement, 2>;
    using cell_ref = SubCell *;
    using cell_ref_list = std::vector<cell_ref>;
    using grid_index = std::tuple<int, int, int>;
    using vec3 = Vec3;
    using level_t = std::uint8_t;
    using particle_index = ParticlesList::particle_index;
    using DataContainer = data::NLDataContainer;

    CellContainer(DataContainer &data, const readdy::model::Context &context,
                  const readdy::util::thread::Config &config);

    virtual ~CellContainer();

    CellContainer(CellContainer &&) = default;

    CellContainer &operator=(CellContainer &&) = default;

    CellContainer(const CellContainer &) = delete;

    CellContainer &operator=(const CellContainer &) = delete;

    sub_cells_t &sub_cells();

    const sub_cells_t &sub_cells() const;

    dimension &n_sub_cells();

    const dimension &n_sub_cells() const;

    const std::size_t n_sub_cells_total() const;

    vec3 &size();

    const vec3 &size() const;

    cell_ref_list &neighbors();

    const cell_ref_list &neighbors() const;

    displacement_arr &maximal_displacements();

    const displacement_arr &maximal_displacements() const;

    cell_index &contiguous_index();

    const cell_index &contiguous_index() const;

    virtual void insert_particle(particle_index index, bool mark_dirty=false) const;

    virtual void insert_particle(particle_index index, bool mark_dirty=false);

    virtual void clear();

    /**
     * return the sub cell that corresponds to the position if it is within bounds (in case of no pbc), otherwise null
     * @todo subtract offset of position so that the index becomes relative
     * @param pos the position
     * @return the corresponding sub cell, otherwise null
     */
    sub_cell *const sub_cell_for_position(const vec3 &pos);

    /**
     * return the sub cell that corresponds to the position if it is within bounds (in case of no pbc), otherwise null
     * @param pos the position
     * @return the corresponding sub cell, otherwise null
     */
    const sub_cell *const sub_cell_for_position(const vec3 &pos) const;

    /**
     * Returns the cell based on the contiguous index, see cell_for_position().
     * @param index the index (triple)
     * @return the cell if found, otherwise null
     */
    sub_cell *const sub_cell_for_grid_index(const grid_index &index);

    const sub_cell *const sub_cell_for_grid_index(const grid_index &index) const;

    sub_cell *const leaf_cell_for_position(const vec3 &pos);

    const sub_cell *const leaf_cell_for_position(const vec3 &pos) const;

    sub_cell *const sub_cell_for_position(const vec3 &pos, level_t level);
    const sub_cell *const sub_cell_for_position(const vec3 &pos, level_t level) const;

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
     *
     * @param max_displacement
     * @return true if everything is alright, false if the displacement of a particle was larger than cutoff + skin
     */
    bool update_sub_cell_displacements_and_mark_dirty(scalar max_cutoff, scalar skin);

    void update_dirty_cells();

    /**
     * sets up subcells with an edge size of >= desired_cell_width
     * @param desired_cell_width
     */
    virtual void subdivide(scalar desired_cell_width);

    virtual void refine_uniformly();

    const DataContainer &data() const;

    DataContainer& data();

    const readdy::model::Context &context() const;

    const readdy::util::thread::Config &config() const;

    const level_t level() const;

    void setup_uniform_neighbors();

    const CellContainer *const root() const;

    CellContainer *const root();

    const vec3& offset() const;

    const CellContainer* const super_cell() const;

    virtual void set_dirty() const;

    virtual void unset_dirty() const;

    virtual void reset_max_displacements();

    void update_root_size();

    /**
     * gives the current number of dirty macro cells, only meaningful after a call to
     * update_sub_cell_displacements_and_mark_dirty().
     * @return the number of dirty macro cells
     */
    const std::size_t n_dirty_macro_cells() const;

    /**
     *
     * @param function
     * @todo this currently only works for leafs at level==2
     */
    virtual void execute_for_each_leaf(const std::function<void(const sub_cell&)> &function);

    void execute_for_each_sub_cell(const std::function<void(sub_cell &)> &function);

    void execute_for_each_sub_cell(const std::function<void(const sub_cell &)> &function) const;

protected:

    CellContainer *_super_cell{nullptr};

    template<typename T, typename Dims>
    static T get_contiguous_index(T i, T j, T k, Dims I, Dims J) {
        return k + j * J + i * I * J;
    }

    sub_cells_t _sub_cells{};
    dimension _n_sub_cells{{0, 0, 0}};
    vec3 _sub_cell_size{0, 0, 0};
    vec3 _size{0, 0, 0};
    vec3 _root_size {0, 0, 0};

    cell_ref_list _neighbors{};
    displacement_arr _maximal_displacements{{0, 0}};
    cell_index _contiguous_index{0};

    // the offset of this cell (lower left corner) w.r.t. the simulation box
    vec3 _offset{0, 0, 0};

    level_t _level{0};

    std::size_t _n_dirty_macro_cells {0};

    DataContainer &_data;
    const readdy::model::Context &_context;
    const readdy::util::thread::Config& _config;
};

}
}
}
}