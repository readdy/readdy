/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file CLLNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include "NeighborList.h"
#include "../model/CPUParticleData.h"
#include "CellLinkedList.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace nl {

class ContiguousCLLNeighborList : public NeighborList {
public:
    ContiguousCLLNeighborList(model::CPUParticleData &data, const readdy::model::KernelContext &context,
                                  const readdy::util::thread::Config &config);

    void set_up(const util::PerformanceNode &node) override;

    void fill_verlet_list(const util::PerformanceNode &node);

    void update(const util::PerformanceNode &node) override;

    void clear(const util::PerformanceNode &node) override;

    void updateData(data_t::update_t &&update) override;

private:
    std::uint8_t cll_radius {1};

    ContiguousCellLinkedList ccll;
    bool _is_set_up {false};
};

}
}
}
}
