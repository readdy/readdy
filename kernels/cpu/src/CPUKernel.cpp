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
 * @file CPUKernel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/cpu/actions/CPUActionFactory.h>
#include <readdy/kernel/cpu/actions/topologies/CPUTopologyActionFactory.h>
#include <readdy/kernel/cpu/observables/CPUObservableFactory.h>

namespace readdy {
namespace kernel {
namespace cpu {
const std::string CPUKernel::name = "CPU";

struct CPUKernel::Impl {
    std::unique_ptr<actions::CPUActionFactory> actionFactory;
    std::unique_ptr<readdy::model::potentials::PotentialFactory> potentialFactory;
    std::unique_ptr<readdy::model::reactions::ReactionFactory> reactionFactory;
    std::unique_ptr<observables::CPUObservableFactory> observableFactory;
    std::unique_ptr<readdy::model::compartments::CompartmentFactory> compartmentFactory;
    std::unique_ptr<CPUStateModel> stateModel;
    std::unique_ptr<readdy::model::KernelContext> context;
    std::unique_ptr<readdy::util::thread::Config> config;
    std::unique_ptr<readdy::model::top::TopologyActionFactory> topologyActionFactory;
};

readdy::model::Kernel* CPUKernel::create() {
    return new CPUKernel();
}

CPUKernel::CPUKernel() : readdy::model::Kernel(name), pimpl(std::make_unique<Impl>()) {
    pimpl->config = std::make_unique<readdy::util::thread::Config>();
    pimpl->config->setMode(readdy::util::thread::ThreadMode::pool);

    pimpl->reactionFactory = std::make_unique<readdy::model::reactions::ReactionFactory>();
    pimpl->context = std::make_unique<readdy::model::KernelContext>();
    pimpl->actionFactory = std::make_unique<actions::CPUActionFactory>(this);
    pimpl->topologyActionFactory = std::make_unique<readdy::kernel::cpu::actions::top::CPUTopologyActionFactory>(this);
    pimpl->stateModel = std::make_unique<CPUStateModel>(pimpl->context.get(), pimpl->config.get(), pimpl->topologyActionFactory.get());
    pimpl->potentialFactory = std::make_unique<readdy::model::potentials::PotentialFactory>();
    pimpl->observableFactory = std::make_unique<observables::CPUObservableFactory>(this);
    pimpl->compartmentFactory = std::make_unique<readdy::model::compartments::CompartmentFactory>();
}

CPUStateModel &CPUKernel::getKernelStateModelInternal() const {
    return *pimpl->stateModel;
}

readdy::model::KernelContext &CPUKernel::getKernelContextInternal() const {
    return *pimpl->context;
}

readdy::model::potentials::PotentialFactory &CPUKernel::getPotentialFactoryInternal() const {
    return *pimpl->potentialFactory;
}

readdy::model::reactions::ReactionFactory &CPUKernel::getReactionFactoryInternal() const {
    return *pimpl->reactionFactory;
}

readdy::model::observables::ObservableFactory &CPUKernel::getObservableFactoryInternal() const {
    return *pimpl->observableFactory;
}

readdy::model::compartments::CompartmentFactory &CPUKernel::getCompartmentFactoryInternal() const {
    return *pimpl->compartmentFactory;
}

unsigned long CPUKernel::getNThreads() const {
    return pimpl->config->nThreads();
}

void CPUKernel::setNThreads(readdy::util::thread::Config::n_threads_type n) {
    pimpl->config->setNThreads(n);
}

readdy::model::actions::ActionFactory &CPUKernel::getActionFactoryInternal() const {
    return *pimpl->actionFactory;
}

readdy::model::top::TopologyActionFactory *CPUKernel::getTopologyActionFactoryInternal() const {
    return pimpl->topologyActionFactory.get();
}

const CPUStateModel &CPUKernel::getCPUKernelStateModel() const {
    return getKernelStateModelInternal();
}

CPUStateModel &CPUKernel::getCPUKernelStateModel() {
    return getKernelStateModelInternal();
}

const readdy::util::thread::Config &CPUKernel::threadConfig() const {
    return *pimpl->config;
}

readdy::util::thread::Config &CPUKernel::threadConfig() {
    return *pimpl->config;
}

void CPUKernel::initialize() {
    readdy::model::Kernel::initialize();
    readdy::conf::cpu::Configuration configuration {};
    const auto &fullConfiguration = getKernelContext().kernelConfiguration();

    if(fullConfiguration.find("CPU") != fullConfiguration.end()) {
        configuration = fullConfiguration["CPU"];
    }

    {
        // thread config
        if (configuration.threadConfig.nThreads > 0) {
            threadConfig().setNThreads(static_cast<unsigned int>(configuration.threadConfig.nThreads));
        }
        threadConfig().setMode(configuration.threadConfig.threadMode);
    }
    {
        // state model config
        getCPUKernelStateModel().configure(configuration);
    }
    for(auto& top : getCPUKernelStateModel().topologies()) {
        top->configure();
        top->updateReactionRates(getKernelContext().topology_registry().structural_reactions_of(top->type()));
    }
}

void CPUKernel::finalize() {
    readdy::model::Kernel::finalize();
    threadConfig().setMode(readdy::util::thread::ThreadMode::inactive);
}

const readdy::util::thread::executor_base &CPUKernel::executor() const {
    return *threadConfig().executor();
}

CPUKernel::~CPUKernel() = default;

}
}
}


const char* name()  {
    return readdy::kernel::cpu::CPUKernel::name.c_str();
}

readdy::model::Kernel* createKernel() {
    return readdy::kernel::cpu::CPUKernel::create();
}
