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
    std::unique_ptr<readdy::model::potentials::PotentialFactory> potentials;
    std::unique_ptr<actions::SCPUActionFactory> actionFactory;
    std::unique_ptr<readdy::model::reactions::ReactionFactory> reactions;
    std::unique_ptr<observables::SCPUObservableFactory> observables;
    std::unique_ptr<model::top::SCPUTopologyActionFactory> topologyActionFactory;
    std::unique_ptr<readdy::model::compartments::CompartmentFactory> compartmentFactory;
};

SCPUKernel::SCPUKernel() : readdy::model::Kernel(name), pimpl(std::make_unique<SCPUKernel::Impl>()) {
    pimpl->actionFactory = std::make_unique<actions::SCPUActionFactory>(this);
    pimpl->topologyActionFactory = std::make_unique<model::top::SCPUTopologyActionFactory>(this);
    pimpl->potentials = std::make_unique<readdy::model::potentials::PotentialFactory>();
    pimpl->reactions = std::make_unique<readdy::model::reactions::ReactionFactory>();
    pimpl->context = std::make_unique<readdy::model::KernelContext>();
    pimpl->model = std::make_unique<SCPUStateModel>(pimpl->context.get(), pimpl->topologyActionFactory.get());
    pimpl->observables = std::make_unique<observables::SCPUObservableFactory>(this);
    pimpl->compartmentFactory = std::make_unique<readdy::model::compartments::CompartmentFactory>();
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

std::vector<std::string> SCPUKernel::getAvailablePotentials() const {
    return pimpl->potentials->getAvailablePotentials();
}

readdy::model::potentials::PotentialFactory &SCPUKernel::getPotentialFactoryInternal() const {
    return *pimpl->potentials;
}

readdy::model::actions::ActionFactory &SCPUKernel::getActionFactoryInternal() const {
    return *pimpl->actionFactory;
}

readdy::model::reactions::ReactionFactory &SCPUKernel::getReactionFactoryInternal() const {
    return *pimpl->reactions;
}

readdy::model::observables::ObservableFactory &SCPUKernel::getObservableFactoryInternal() const {
    return *pimpl->observables;
}

readdy::model::top::TopologyActionFactory *SCPUKernel::getTopologyActionFactoryInternal() const {
    return pimpl->topologyActionFactory.get();
}
readdy::model::compartments::CompartmentFactory &SCPUKernel::getCompartmentFactoryInternal() const {
    return *pimpl->compartmentFactory;
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
        top->updateReactionRates(getKernelContext().topology_registry().reactions_of(top->type()));
    }
}

SCPUKernel &SCPUKernel::operator=(SCPUKernel &&rhs) noexcept = default;

SCPUKernel::SCPUKernel(SCPUKernel &&rhs) noexcept = default;

}
}
}




