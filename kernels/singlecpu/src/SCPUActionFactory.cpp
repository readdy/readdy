/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


//
// Created by clonker on 08.04.16.
//

#include <memory>

#include <readdy/kernel/singlecpu/actions/SCPUActionFactory.h>
#include <readdy/kernel/singlecpu/actions/SCPUEulerBDIntegrator.h>
#include <readdy/kernel/singlecpu/actions/SCPUCalculateForces.h>
#include <readdy/kernel/singlecpu/actions/SCPUReactionImpls.h>
#include <readdy/kernel/singlecpu/actions/SCPUUpdateNeighborList.h>
#include <readdy/kernel/singlecpu/actions/SCPUEvaluateCompartments.h>
#include <readdy/kernel/singlecpu/actions/SCPUEvaluateTopologyReactions.h>

namespace core_actions = readdy::model::actions;

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
SCPUActionFactory::SCPUActionFactory(SCPUKernel *const kernel) : kernel(kernel) {}

namespace rma = readdy::model::actions;

std::vector<std::string> SCPUActionFactory::getAvailableActions() const {
    return {
            rma::getActionName<rma::AddParticles>(), rma::getActionName<rma::EulerBDIntegrator>(),
            rma::getActionName<rma::CalculateForces>(),
            rma::getActionName<rma::UpdateNeighborList>(),
            rma::getActionName<rma::reactions::UncontrolledApproximation>(),
            rma::getActionName<rma::reactions::Gillespie>(),
            rma::getActionName<rma::top::EvaluateTopologyReactions>()
    };
}

std::unique_ptr<readdy::model::actions::EulerBDIntegrator> SCPUActionFactory::eulerBDIntegrator(scalar timeStep) const {
    return {std::make_unique<SCPUEulerBDIntegrator>(kernel, timeStep)};
}

std::unique_ptr<readdy::model::actions::CalculateForces> SCPUActionFactory::calculateForces() const {
    return {std::make_unique<SCPUCalculateForces>(kernel)};
}

std::unique_ptr<readdy::model::actions::UpdateNeighborList>
SCPUActionFactory::updateNeighborList(readdy::model::actions::UpdateNeighborList::Operation operation,
                                      scalar skinSize) const {
    return {std::make_unique<SCPUUpdateNeighborList>(kernel, operation, skinSize)};
}

std::unique_ptr<readdy::model::actions::EvaluateCompartments> SCPUActionFactory::evaluateCompartments() const {
    return {std::make_unique<SCPUEvaluateCompartments>(kernel)};
}

std::unique_ptr<readdy::model::actions::reactions::UncontrolledApproximation>
SCPUActionFactory::uncontrolledApproximation(scalar timeStep) const {
    return {std::make_unique<reactions::SCPUUncontrolledApproximation>(kernel, timeStep)};
}

std::unique_ptr<readdy::model::actions::reactions::Gillespie> SCPUActionFactory::gillespie(scalar timeStep) const {
    return {std::make_unique<reactions::SCPUGillespie>(kernel, timeStep)};
}

std::unique_ptr<readdy::model::actions::top::EvaluateTopologyReactions>
SCPUActionFactory::evaluateTopologyReactions(scalar timeStep) const {
    return {std::make_unique<top::SCPUEvaluateTopologyReactions>(kernel, timeStep)};
}

std::unique_ptr<readdy::model::actions::reactions::DetailedBalance>
SCPUActionFactory::detailedBalance(scalar timeStep) const {
    return {std::make_unique<reactions::SCPUDetailedBalance>(kernel, timeStep)};
}

}
}
}
}


