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
namespace cpu {
namespace nl {

namespace detail {
class DirtyFlag {
public:
    DirtyFlag() = default;
    DirtyFlag(DirtyFlag&&);
    DirtyFlag& operator=(DirtyFlag&&);
    DirtyFlag(const DirtyFlag&) = delete;
    DirtyFlag& operator=(const DirtyFlag&) = delete;

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

    SubCell(CellContainer *const super_cell, const vec3 &offset);

    SubCell(SubCell&& rhs) = default;

    const bool is_leaf() const;

    virtual void update_displacements() override;

    virtual void subdivide(const scalar desired_cell_width) override;

    virtual void reset_max_displacements() override;

    virtual void refine_uniformly() override;

    void setup_uniform_neighbors(const std::uint8_t radius);

    void insert_particle(const particle_index index) const override;

    void insert_particle(const particle_index index);

    virtual void clear() override;

    const ParticlesList& particles() const;

    const bool is_dirty() const;

    void set_dirty() const;

    void unset_dirty() const;

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