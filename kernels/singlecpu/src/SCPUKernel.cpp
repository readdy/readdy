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


//
// Created by clonker on 07.03.16.
//

#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
const std::string SCPUKernel::name = "SingleCPU";

SCPUKernel::SCPUKernel() : readdy::model::Kernel(name), _actionFactory(this), _topologyActionFactory(this),
                           _model(_context, &_topologyActionFactory), _observables(this) {}

/**
 * factory method
 */
std::unique_ptr<SCPUKernel> SCPUKernel::create() {
    return std::make_unique<SCPUKernel>();
}

/**
 * Destructor: default
 */
SCPUKernel::~SCPUKernel() = default;

void SCPUKernel::initialize() {
    readdy::model::Kernel::initialize();
    for(auto& top : getSCPUKernelStateModel().topologies()) {
        top->configure();
        top->updateReactionRates(context().topologyRegistry().structuralReactionsOf(top->type()));
    }
    getSCPUKernelStateModel().reactionRecords().clear();
    getSCPUKernelStateModel().resetReactionCounts();
    getSCPUKernelStateModel().virial() = Matrix33{{{0, 0, 0, 0, 0, 0, 0, 0, 0}}};
}

}
}
}




