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
 * << detailed description >>
 *
 * @file CellDecompositionNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 01.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/kernel/cpu/util/config.h>

#include "NeighborList.h"
#include "CellContainer.h"
#include "SubCell.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class CellDecompositionNeighborList : public NeighborList {
public:

    CellDecompositionNeighborList(data::EntryDataContainer *data, const readdy::model::Context &context,
                                  const readdy::util::thread::Config &config);

    CellDecompositionNeighborList(const readdy::model::Context &context,
                                  const readdy::util::thread::Config &config);

    bool is_adaptive() const override;

    void set_up(const util::PerformanceNode &node) override;

    void update(const util::PerformanceNode &node) override;

    void clear(const util::PerformanceNode &node) override;

    void updateData(DataUpdate &&update) override;

    void fill_container();

    std::size_t size() const override;

    void fill_verlet_list();

    void fill_cell_verlet_list(const CellContainer::sub_cell &sub_cell);

    virtual const_iterator cbegin() const override;

    virtual const_iterator cend() const override;

    const data::EntryDataContainer *data() const override;

    data::EntryDataContainer *data() override;

private:
    bool _is_set_up {false};

    data::NLDataContainer _data;
    CellContainer _cell_container;
};

}
}
}
}
