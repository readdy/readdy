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
 * @file SingleCPUUpdateNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.07.16
 */
#include <readdy/kernel/singlecpu/actions/SCPUUpdateNeighborList.h>

namespace core_actions = readdy::model::actions;

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {

void SCPUUpdateNeighborList::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    switch (operation) {
        case init:
            kernel->getSCPUKernelStateModel().getNeighborList()->setUp(skinSize > 0 ? skinSize : 0, 1, node);
            break;
        case clear:
            kernel->stateModel().clearNeighborList();
            break;
        case update:
            kernel->getSCPUKernelStateModel().getNeighborList()->update(node);
            break;
    }
}

SCPUUpdateNeighborList::SCPUUpdateNeighborList(SCPUKernel *const kernel, core_actions::UpdateNeighborList::Operation op,
                                               scalar skinSize)
        : UpdateNeighborList(op, skinSize), kernel(kernel){
}

}
}
}
}
