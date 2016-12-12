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
#include <readdy/kernel/singlecpu/programs/SCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SCPUTestProgram.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialFactory.h>
#include <readdy/kernel/singlecpu/reactions/SCPUReactionFactory.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservableFactory.h>


namespace readdy {
namespace kernel {
namespace scpu {
const std::string SCPUKernel::name = "SingleCPU";
struct SCPUKernel::Impl {
    std::unique_ptr<readdy::model::KernelContext> context;
    std::unique_ptr<SCPUStateModel> model;
    std::unique_ptr<potentials::SCPUPotentialFactory> potentials;
    std::unique_ptr<programs::SCPUProgramFactory> programs;
    std::unique_ptr<reactions::SCPUReactionFactory> reactions;
    std::unique_ptr<observables::SCPUObservableFactory> observables;
};

SCPUKernel::SCPUKernel() : readdy::model::Kernel(name), pimpl(std::make_unique<SCPUKernel::Impl>()) {
    pimpl->programs = std::make_unique<programs::SCPUProgramFactory>(this);
    pimpl->potentials = std::make_unique<potentials::SCPUPotentialFactory>(this);
    pimpl->reactions = std::make_unique<reactions::SCPUReactionFactory>(this);
    pimpl->context = std::make_unique<readdy::model::KernelContext>();
    pimpl->model = std::make_unique<SCPUStateModel>(pimpl->context.get());
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

SCPUStateModel &SCPUKernel::getKernelStateModel() const {
    return *pimpl->model;
}

readdy::model::KernelContext &SCPUKernel::getKernelContext() const {
    return *pimpl->context;
}

std::vector<std::string> SCPUKernel::getAvailablePotentials() const {
    return pimpl->potentials->getAvailablePotentials();
}

std::unique_ptr<readdy::model::potentials::Potential> SCPUKernel::createPotential(std::string &name) const {
    return pimpl->potentials->createPotential(name);
}

readdy::model::potentials::PotentialFactory &SCPUKernel::getPotentialFactory() const {
    return *pimpl->potentials;
}

readdy::model::programs::ProgramFactory &SCPUKernel::getProgramFactory() const {
    return *pimpl->programs;
}

readdy::model::reactions::ReactionFactory &SCPUKernel::getReactionFactory() const {
    return *pimpl->reactions;
}

readdy::model::observables::ObservableFactory &SCPUKernel::getObservableFactory() const {
    return *pimpl->observables;
}


SCPUKernel &SCPUKernel::operator=(SCPUKernel &&rhs) = default;

SCPUKernel::SCPUKernel(SCPUKernel &&rhs) = default;

}
}
}



