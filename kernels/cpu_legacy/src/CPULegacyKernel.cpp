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
 * @file CPULegacyKernel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu_legacy/CPULegacyKernel.h>
#include <readdy/kernel/cpu_legacy/actions/CPUActionFactory.h>
#include <readdy/kernel/cpu_legacy/actions/topologies/CPUTopologyActionFactory.h>
#include <readdy/kernel/cpu_legacy/observables/CPUObservableFactory.h>

namespace readdy {
namespace kernel {
namespace cpu_legacy {
const std::string CPULegacyKernel::name = "CPU_Legacy";

struct CPULegacyKernel::Impl {
    std::unique_ptr<actions::CPUActionFactory> actionFactory;
    std::unique_ptr<observables::CPUObservableFactory> observableFactory;
    std::unique_ptr<CPUStateModel> stateModel;
    std::unique_ptr<readdy::util::thread::Config> config;
    std::unique_ptr<readdy::model::top::TopologyActionFactory> topologyActionFactory;
};

readdy::model::Kernel *CPULegacyKernel::create() {
    return new CPULegacyKernel();
}

CPULegacyKernel::CPULegacyKernel() : readdy::model::Kernel(name), pimpl(std::make_unique<Impl>()) {
    pimpl->config = std::make_unique<readdy::util::thread::Config>();
    pimpl->config->setMode(readdy::util::thread::ThreadMode::pool);

    pimpl->actionFactory = std::make_unique<actions::CPUActionFactory>(this);
    pimpl->topologyActionFactory = std::make_unique<actions::top::CPUTopologyActionFactory>(this);
    pimpl->stateModel = std::make_unique<CPUStateModel>(_context, pimpl->config.get(),
                                                        pimpl->topologyActionFactory.get());
    pimpl->observableFactory = std::make_unique<observables::CPUObservableFactory>(this);
}

CPUStateModel &CPULegacyKernel::getKernelStateModelInternal() const {
    return *pimpl->stateModel;
}

unsigned long CPULegacyKernel::getNThreads() const {
    return pimpl->config->nThreads();
}

void CPULegacyKernel::setNThreads(readdy::util::thread::Config::n_threads_type n) {
    pimpl->config->setNThreads(n);
}

const CPUStateModel &CPULegacyKernel::getCPULegacyKernelStateModel() const {
    return getKernelStateModelInternal();
}

CPUStateModel &CPULegacyKernel::getCPULegacyKernelStateModel() {
    return getKernelStateModelInternal();
}

const readdy::util::thread::Config &CPULegacyKernel::threadConfig() const {
    return *pimpl->config;
}

readdy::util::thread::Config &CPULegacyKernel::threadConfig() {
    return *pimpl->config;
}

void CPULegacyKernel::initialize() {
    readdy::model::Kernel::initialize();

    const auto &fullConfiguration = context().kernelConfiguration();

    const auto &configuration = fullConfiguration.cpuLegacy;
    {
        // thread config
        if (configuration.threadConfig.nThreads > 0) {
            threadConfig().setNThreads(static_cast<unsigned int>(configuration.threadConfig.nThreads));
        }
        threadConfig().setMode(configuration.threadConfig.threadMode);
    }
    {
        // state model config
        getCPULegacyKernelStateModel().configure(configuration);
    }
    for (auto &top : getCPULegacyKernelStateModel().topologies()) {
        top->configure();
        top->updateReactionRates(context().topology_registry().structuralReactionsOf(top->type()));
    }
    getCPULegacyKernelStateModel().reactionRecords().clear();
    getCPULegacyKernelStateModel().resetReactionCounts();
}

void CPULegacyKernel::finalize() {
    readdy::model::Kernel::finalize();
    threadConfig().setMode(readdy::util::thread::ThreadMode::inactive);
}

const readdy::util::thread::executor_base &CPULegacyKernel::executor() const {
    return *threadConfig().executor();
}

const model::StateModel &CPULegacyKernel::stateModel() const {
    return *pimpl->stateModel;
}

model::StateModel &CPULegacyKernel::stateModel() {
    return *pimpl->stateModel;
}

const model::actions::ActionFactory &CPULegacyKernel::getActionFactory() const {
    return *pimpl->actionFactory;
}

model::actions::ActionFactory &CPULegacyKernel::getActionFactory() {
    return *pimpl->actionFactory;
}

const model::top::TopologyActionFactory *const CPULegacyKernel::getTopologyActionFactory() const {
    return pimpl->topologyActionFactory.get();
}

model::top::TopologyActionFactory *const CPULegacyKernel::getTopologyActionFactory() {
    return pimpl->topologyActionFactory.get();
}

const model::observables::ObservableFactory &CPULegacyKernel::getObservableFactory() const {
    return *pimpl->observableFactory;
}

model::observables::ObservableFactory &CPULegacyKernel::getObservableFactory() {
    return *pimpl->observableFactory;
}

CPULegacyKernel::~CPULegacyKernel() = default;

}
}
}


const char *name() {
    return readdy::kernel::cpu_legacy::CPULegacyKernel::name.c_str();
}

readdy::model::Kernel *createKernel() {
    return readdy::kernel::cpu_legacy::CPULegacyKernel::create();
}
