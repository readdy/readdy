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

SCPUKernel::SCPUKernel() : readdy::model::Kernel(name) {
    _actionFactory = std::make_unique<actions::SCPUActionFactory>(this);
    _topologyActionFactory = std::make_unique<model::top::SCPUTopologyActionFactory>(this);
    _model = std::make_unique<SCPUStateModel>(_context, _topologyActionFactory.get());
    _observables = std::make_unique<observables::SCPUObservableFactory>(this);
}

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

SCPUStateModel &SCPUKernel::getKernelStateModelInternal() const {
    return *_model;
}

readdy::model::actions::ActionFactory &SCPUKernel::getActionFactoryInternal() const {
    return *_actionFactory;
}

readdy::model::observables::ObservableFactory &SCPUKernel::getObservableFactoryInternal() const {
    return *_observables;
}

readdy::model::top::TopologyActionFactory *SCPUKernel::getTopologyActionFactoryInternal() const {
    return _topologyActionFactory.get();
}

const SCPUStateModel &SCPUKernel::getSCPUKernelStateModel() const {
    return getKernelStateModelInternal();
}

SCPUStateModel &SCPUKernel::getSCPUKernelStateModel() {
    return getKernelStateModelInternal();
}

void SCPUKernel::initialize() {
    readdy::model::Kernel::initialize();
    for(auto& top : getSCPUKernelStateModel().topologies()) {
        top->configure();
        top->updateReactionRates(context().topology_registry().structuralReactionsOf(top->type()));
    }
    getSCPUKernelStateModel().resetReactionCounts();
}

}
}
}




