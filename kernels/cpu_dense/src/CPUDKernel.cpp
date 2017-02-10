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
 * @file Kernel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#include "readdy/kernel/cpu_dense/CPUDKernel.h"
#include <readdy/kernel/cpu_dense/observables/CPUDObservableFactory.h>
#include <readdy/kernel/cpu_dense/actions/CPUDActionFactory.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {

const std::string CPUDKernel::name = "CPU_Dense";

struct CPUDKernel::Impl {
    std::unique_ptr<actions::CPUDActionFactory> actionFactory;
    std::unique_ptr<readdy::model::potentials::PotentialFactory> potentialFactory;
    std::unique_ptr<readdy::model::reactions::ReactionFactory> reactionFactory;
    std::unique_ptr<observables::CPUDObservableFactory> observableFactory;
    std::unique_ptr<readdy::model::compartments::CompartmentFactory> compartmentFactory;
    std::unique_ptr<CPUDStateModel> stateModel;
    std::unique_ptr<readdy::model::KernelContext> context;
    std::unique_ptr<readdy::util::thread::Config> config;
};

readdy::model::Kernel* CPUDKernel::create() {
    return new CPUDKernel();
}


CPUDKernel::CPUDKernel() : readdy::model::Kernel(name), pimpl(std::make_unique<Impl>()) {
    pimpl->config = std::make_unique<readdy::util::thread::Config>();
    pimpl->reactionFactory = std::make_unique<readdy::model::reactions::ReactionFactory>();
    pimpl->context = std::make_unique<readdy::model::KernelContext>();
    pimpl->actionFactory = std::make_unique<actions::CPUDActionFactory>(this);
    pimpl->stateModel = std::make_unique<CPUDStateModel>(pimpl->context.get(), pimpl->config.get());
    pimpl->potentialFactory = std::make_unique<readdy::model::potentials::PotentialFactory>();
    pimpl->observableFactory = std::make_unique<observables::CPUDObservableFactory>(this);
    pimpl->compartmentFactory = std::make_unique<readdy::model::compartments::CompartmentFactory>();
}

readdy::model::actions::ActionFactory &CPUDKernel::getActionFactoryInternal() const {
    return *pimpl->actionFactory;
}


CPUDStateModel &CPUDKernel::getKernelStateModelInternal() const {
    return *pimpl->stateModel;
}

readdy::model::KernelContext &CPUDKernel::getKernelContextInternal() const {
    return *pimpl->context;
}

readdy::model::potentials::PotentialFactory &CPUDKernel::getPotentialFactoryInternal() const {
    return *pimpl->potentialFactory;
}

readdy::model::reactions::ReactionFactory &CPUDKernel::getReactionFactoryInternal() const {
    return *pimpl->reactionFactory;
}

readdy::model::observables::ObservableFactory &CPUDKernel::getObservableFactoryInternal() const {
    return *pimpl->observableFactory;
}

readdy::model::compartments::CompartmentFactory &CPUDKernel::getCompartmentFactoryInternal() const {
    return *pimpl->compartmentFactory;
}

unsigned long CPUDKernel::getNThreads() const {
    return pimpl->config->nThreads();
}

void CPUDKernel::setNThreads(readdy::util::thread::Config::n_threads_t n) {
    pimpl->config->setNThreads(n);
}

readdy::model::top::TopologyActionFactory *CPUDKernel::getTopologyActionFactoryInternal() const {
    return nullptr;
}

const CPUDStateModel &CPUDKernel::getCPUDKernelStateModel() const {
    return getKernelStateModelInternal();
}

CPUDStateModel &CPUDKernel::getCPUDKernelStateModel() {
    return getKernelStateModelInternal();
}

CPUDKernel::~CPUDKernel() = default;

}
}
}


const char* name()  {
    return readdy::kernel::cpu_dense::CPUDKernel::name.c_str();
}

readdy::model::Kernel* createKernel() {
    return readdy::kernel::cpu_dense::CPUDKernel::create();
}
