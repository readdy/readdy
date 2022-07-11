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
 * @file CPUActionFactory.h
 * @brief Declaration of CPU action factory
 * @author clonker
 * @date 23.06.16
 */

#pragma once

#include <readdy/model/actions/ActionFactory.h>

namespace readdy::kernel::cpu {

class CPUKernel;
namespace actions {
class CPUActionFactory : public readdy::model::actions::ActionFactory {
    CPUKernel *const kernel;
public:
    explicit CPUActionFactory(CPUKernel *kernel);

    std::unique_ptr<model::actions::AddParticles>
    addParticles(const std::vector<model::Particle> &particles) const override;

    std::unique_ptr<model::actions::EulerBDIntegrator> eulerBDIntegrator(scalar timeStep) const override;

    std::unique_ptr<readdy::model::actions::CalculateForces> calculateForces() const override;

    std::unique_ptr<model::actions::CreateNeighborList> createNeighborList(scalar interactionDistance) const override;

    std::unique_ptr<model::actions::UpdateNeighborList> updateNeighborList() const override;

    std::unique_ptr<model::actions::ClearNeighborList> clearNeighborList() const override;

    std::unique_ptr<model::actions::EvaluateCompartments> evaluateCompartments() const override;

    std::unique_ptr<model::actions::reactions::UncontrolledApproximation>
    uncontrolledApproximation(scalar timeStep) const override;

    std::unique_ptr<model::actions::reactions::Gillespie>
    gillespie(scalar timeStep) const override;

    std::unique_ptr<model::actions::reactions::DetailedBalance>
    detailedBalance(scalar timeStep) const override;

    std::unique_ptr<model::actions::top::EvaluateTopologyReactions>
    evaluateTopologyReactions(scalar timeStep) const override;

    std::unique_ptr<model::actions::top::BreakBonds>
    breakBonds(scalar timeStep, readdy::model::actions::top::BreakConfig config) const override;

    std::unique_ptr<model::actions::EvaluateObservables> evaluateObservables() const override;

    std::unique_ptr<model::actions::MakeCheckpoint>
    makeCheckpoint(std::string base, std::size_t maxNSaves) const override;

    std::unique_ptr<model::actions::InitializeKernel> initializeKernel() const override;
};

}
}
