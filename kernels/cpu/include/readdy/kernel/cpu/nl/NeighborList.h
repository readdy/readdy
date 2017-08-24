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
 * @file NeighborList.h
 * @brief 
 * @author clonker
 * @date 4/21/17
 */
#pragma once

#include <readdy/common/macros.h>
#include <readdy/kernel/cpu/util/config.h>

#include "CellContainer.h"
#include "SubCell.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class NeighborList {
public:
    using skin_size_t = scalar;
    using data_t = readdy::kernel::cpu::model::CPUParticleData;
    using neighbors_t = data_t::neighbors_t;
    using iterator = data_t::neighbors_list_iterator;
    using const_iterator = data_t::neighbors_list_const_iterator;

    NeighborList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                 const readdy::util::thread::Config &config, bool adaptive=true, skin_size_t skin = 0,
                 bool hilbert_sort = true);

    void set_up();

    void update();

    void clear();

    void clear_cells();

    const model::CPUParticleData& data() const;

    const readdy::util::thread::Config& config() const;

    const CellContainer& cell_container() const;

    bool& adaptive();

    const bool& adaptive() const;

    bool& performs_hilbert_sort();

    const bool& performs_hilbert_sort() const;

    void sort_by_hilbert_curve();

    void updateData(data_t::update_t &&update);

    void displace(data_t::iterator iter, const readdy::model::Vec3 &vec);

    void displace(data_t::Entry &entry, const readdy::model::Vec3 &delta);

    void displace(data_t::index_t entry, const readdy::model::Vec3 &delta);

    iterator begin();

    iterator end();

    const_iterator cbegin() const;

    const_iterator cend() const;

    const neighbors_t &neighbors_of(const data_t::index_t entry) const;

    skin_size_t &skin();

    const skin_size_t &skin() const;

private:

    void fill_container();

    /**
     * should be called once containers are all valid / filled
     */
    void fill_verlet_list();

    void fill_cell_verlet_list(const CellContainer::sub_cell &sub_cell, bool reset_displacement);

    void handle_dirty_cells();

    CellContainer _cell_container;

    skin_size_t _skin;
    scalar _max_cutoff;
    scalar _max_cutoff_skin_squared {0};
    bool _hilbert_sort {true};
    bool _adaptive {true};
    bool _is_set_up {false};
    model::CPUParticleData &_data;
    const readdy::model::KernelContext &_context;
    const readdy::util::thread::Config &_config;
};

}
}
}
}