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
#include <readdy/kernel/singlecpu/programs/SCPUUpdateNeighborList.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace programs {
SCPUUpdateNeighborList::SCPUUpdateNeighborList(SCPUKernel *kernel) : kernel(kernel) {}

void SCPUUpdateNeighborList::execute() {
    switch (action) {
        case create:
            kernel->getKernelStateModel().updateNeighborList();
            break;
        case clear:
            kernel->getKernelStateModel().clearNeighborList();
            break;
    }
}
}
}
}
}
