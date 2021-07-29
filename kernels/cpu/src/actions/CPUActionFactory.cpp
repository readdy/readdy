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


/**
 * @file CPUActionFactory.cpp
 * @brief CPU kernel implementation of Action factory
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/actions/CPUActionFactory.h>
#include <readdy/kernel/cpu/actions/CPUEulerBDIntegrator.h>
#include <readdy/kernel/cpu/actions/CPUCreateNeighborList.h>
#include <readdy/kernel/cpu/actions/CPUCalculateForces.h>
#include <readdy/kernel/cpu/actions/CPUEvaluateCompartments.h>
#include <readdy/kernel/cpu/actions/reactions/CPUGillespie.h>
#include <readdy/kernel/cpu/actions/reactions/CPUUncontrolledApproximation.h>
#include <readdy/kernel/cpu/actions/CPUEvaluateTopologyReactions.h>
#include <readdy/kernel/cpu/actions/CPUBreakBonds.h>
#include <readdy/kernel/cpu/actions/CPUMiscActions.h>

namespace core_p = readdy::model::actions;

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
CPUActionFactory::CPUActionFactory(CPUKernel *const kernel) : kernel(kernel) { }

std::unique_ptr<model::actions::AddParticles>
CPUActionFactory::addParticles(const std::vector<model::Particle> &particles) const {
    return {std::make_unique<readdy::model::actions::AddParticles>(kernel, particles)};
}

std::unique_ptr<model::actions::EulerBDIntegrator> CPUActionFactory::eulerBDIntegrator(scalar timeStep) const {
    return {std::make_unique<CPUEulerBDIntegrator>(kernel, timeStep)};
}

std::unique_ptr<readdy::model::actions::CalculateForces> CPUActionFactory::calculateForces() const {
    return {std::make_unique<CPUCalculateForces>(kernel)};
}

std::unique_ptr<model::actions::CreateNeighborList>
CPUActionFactory::createNeighborList(scalar cutoffDistance) const {
    return {std::make_unique<CPUCreateNeighborList>(kernel, cutoffDistance)};
}

std::unique_ptr<model::actions::UpdateNeighborList> CPUActionFactory::updateNeighborList() const {
    return {std::make_unique<CPUUpdateNeighborList>(kernel)};
}

std::unique_ptr<model::actions::ClearNeighborList> CPUActionFactory::clearNeighborList() const {
    return {std::make_unique<CPUClearNeighborList>(kernel)};
}

std::unique_ptr<model::actions::EvaluateCompartments> CPUActionFactory::evaluateCompartments() const {
    return {std::make_unique<CPUEvaluateCompartments>(kernel)};
}

std::unique_ptr<model::actions::reactions::UncontrolledApproximation>
CPUActionFactory::uncontrolledApproximation(scalar timeStep) const {
    return {std::make_unique<reactions::CPUUncontrolledApproximation>(kernel, timeStep)};
}

std::unique_ptr<model::actions::reactions::Gillespie>
CPUActionFactory::gillespie(scalar timeStep) const {
    return {std::make_unique<reactions::CPUGillespie>(kernel, timeStep)};
}

std::unique_ptr<model::actions::top::EvaluateTopologyReactions>
CPUActionFactory::evaluateTopologyReactions(scalar timeStep) const {
    return {std::make_unique<top::CPUEvaluateTopologyReactions>(kernel, timeStep)};
}

std::unique_ptr<model::actions::reactions::DetailedBalance>
CPUActionFactory::detailedBalance(scalar timeStep) const {
    throw std::invalid_argument("DetailedBalance reaction handler not implemented for CPU");
}

std::unique_ptr<model::actions::top::BreakBonds>
CPUActionFactory::breakBonds(scalar timeStep, readdy::model::actions::top::BreakConfig config) const {
    return {std::make_unique<top::CPUBreakBonds>(kernel, timeStep, config)};
}

std::unique_ptr<model::actions::EvaluateObservables> CPUActionFactory::evaluateObservables() const {
    return {std::make_unique<CPUEvaluateObservables>(kernel)};
}

std::unique_ptr<model::actions::MakeCheckpoint>
CPUActionFactory::makeCheckpoint(std::string base, std::size_t maxNSaves, std::string checkpointFormat) const {
    return {std::make_unique<CPUMakeCheckpoint>(kernel, base, maxNSaves, checkpointFormat)};
}

std::unique_ptr<model::actions::InitializeKernel> CPUActionFactory::initializeKernel() const {
    return {std::make_unique<CPUInitializeKernel>(kernel)};
}

}
}
}
}
