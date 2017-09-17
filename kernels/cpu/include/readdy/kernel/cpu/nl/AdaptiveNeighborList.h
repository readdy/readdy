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

#include <readdy/kernel/cpu/util/config.h>

#include <readdy/kernel/cpu/data/NLDataContainer.h>

#include "CellContainer.h"
#include "SubCell.h"
#include "NeighborList.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class AdaptiveNeighborList : public NeighborList{
public:
    using skin_size_t = scalar;
    using data_type = data::NLDataContainer;

    AdaptiveNeighborList(data::EntryDataContainer *data, const readdy::model::KernelContext &context, const readdy::util::thread::Config &config,
                         bool hilbert_sort = true);

    AdaptiveNeighborList(const readdy::model::KernelContext &context, const readdy::util::thread::Config &config,
                         bool hilbert_sort = true);

    ~AdaptiveNeighborList() = default;

    bool is_adaptive() const override;

    void set_up(const util::PerformanceNode &node) override;

    void update(const util::PerformanceNode &node) override;

    void clear(const util::PerformanceNode &node) override;

    void clear_cells();

    const readdy::util::thread::Config& config() const;

    const CellContainer& cell_container() const;

    bool& performs_hilbert_sort();

    const bool& performs_hilbert_sort() const;

    void sort_by_hilbert_curve();

    void updateData(data_type::DataUpdate &&update) override;

    virtual const neighbors_type &neighbors_of(const data_type::size_type entry) const override;

    virtual const_iterator cbegin() const override;

    virtual const_iterator cend() const override;

    data::EntryDataContainer *data() override;

    const data::EntryDataContainer * data() const override;

private:

    void fill_container();

    /**
     * should be called once containers are all valid / filled
     */
    void fill_verlet_list();

    void fill_cell_verlet_list(const CellContainer::sub_cell &sub_cell, bool reset_displacement);

    void handle_dirty_cells();

    bool _hilbert_sort {true};
    bool _is_set_up {false};

    data::NLDataContainer _data;

    CellContainer _cell_container;
};

}
}
}
}