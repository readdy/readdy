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
#include <readdy/kernel/singlecpu/actions/SCPUActionFactory.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservableFactory.h>
#include <readdy/kernel/singlecpu/model/topologies/SCPUTopologyActionFactory.h>


namespace readdy {
namespace kernel {
namespace scpu {
const std::string SCPUKernel::name = "SingleCPU";
struct SCPUKernel::Impl {
    std::unique_ptr<readdy::model::KernelContext> context;
    std::unique_ptr<SCPUStateModel> model;
    std::unique_ptr<actions::SCPUActionFactory> actionFactory;
    std::unique_ptr<observables::SCPUObservableFactory> observables;
    std::unique_ptr<model::top::SCPUTopologyActionFactory> topologyActionFactory;
};

SCPUKernel::SCPUKernel() : readdy::model::Kernel(name), pimpl(std::make_unique<SCPUKernel::Impl>()) {
    pimpl->actionFactory = std::make_unique<actions::SCPUActionFactory>(this);
    pimpl->topologyActionFactory = std::make_unique<model::top::SCPUTopologyActionFactory>(this);
    pimpl->context = std::make_unique<readdy::model::KernelContext>();
    pimpl->model = std::make_unique<SCPUStateModel>(pimpl->context.get(), pimpl->topologyActionFactory.get());
    pimpl->observables = std::make_unique<observables::SCPUObservableFactory>(this);
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
    return *pimpl->model;
}

readdy::model::KernelContext &SCPUKernel::getKernelContextInternal() const {
    return *pimpl->context;
}

readdy::model::actions::ActionFactory &SCPUKernel::getActionFactoryInternal() const {
    return *pimpl->actionFactory;
}

readdy::model::observables::ObservableFactory &SCPUKernel::getObservableFactoryInternal() const {
    return *pimpl->observables;
}

readdy::model::top::TopologyActionFactory *SCPUKernel::getTopologyActionFactoryInternal() const {
    return pimpl->topologyActionFactory.get();
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
        top->updateReactionRates(getKernelContext().topology_registry().structural_reactions_of(top->type()));
    }
}

SCPUKernel &SCPUKernel::operator=(SCPUKernel &&rhs) noexcept = default;

SCPUKernel::SCPUKernel(SCPUKernel &&rhs) noexcept = default;

}
}
}




