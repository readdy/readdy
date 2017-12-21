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

#include <atomic>

#include "CellContainer.h"

namespace readdy {
namespace kernel {
namespace cpu_legacy {
namespace nl {

namespace detail {
class DirtyFlag {
public:
    DirtyFlag() = default;
    DirtyFlag(DirtyFlag&&) noexcept;
    DirtyFlag& operator=(DirtyFlag&&) noexcept;
    DirtyFlag(const DirtyFlag&) = delete;
    DirtyFlag& operator=(const DirtyFlag&) = delete;
    ~DirtyFlag() = default;

    void set() const;

    void unset() const;

    bool get() const;

private:
    mutable std::atomic<bool> _is_dirty{false};
};
}

class SubCell : public CellContainer {
    using super = CellContainer;
public:
    using particle_ref = int;

    SubCell(CellContainer* super_cell, const vec3 &offset);

    SubCell(SubCell&& rhs) = default;

    SubCell& operator=(SubCell&&) = default;

    SubCell(const SubCell&) = delete;

    SubCell& operator=(const SubCell&) = delete;

    ~SubCell() override = default;

    const bool is_leaf() const;

    void update_displacements() override;

    void subdivide(scalar desired_cell_width) override;

    void reset_max_displacements() override;

    void refine_uniformly() override;

    void setup_uniform_neighbors(std::uint8_t radius);

    void insert_particle(particle_index index, bool mark_dirty=false) const override;

    void insert_particle(particle_index index, bool mark_dirty=false) override;

    //virtual void execute_for_each_leaf(const std::function<void(const SubCell&)> &function) override;

    void clear() override;

    const ParticlesList& particles() const;

    const bool is_dirty() const;

    void set_dirty() const override;

    void unset_dirty() const override;

    const bool neighbor_dirty() const;

    void reset_particles_displacements();

    ParticlesList::particle_indices collect_contained_particles() const;

private:
    bool _is_leaf{true};

    detail::DirtyFlag _dirty_flag {};
    ParticlesList _particles_list {};

    // change visibility to private, this should not be used with sub cells
    void update_sub_cell_displacements() override;

};

}
}
}
}