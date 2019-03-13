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
 * @file SCPUActionFactory.h
 * @brief Declaration of the Single CPU action factory.
 * @author clonker
 * @author chrisfroe
 * @date 08.04.16
 */

#pragma once

#include <readdy/model/actions/ActionFactory.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>

namespace readdy {
namespace kernel {
namespace scpu {
class SCPUKernel;
namespace actions {
class SCPUActionFactory : public readdy::model::actions::ActionFactory {
public:
    explicit SCPUActionFactory(SCPUKernel* kernel);

    std::vector<std::string> getAvailableActions() const override;

    std::unique_ptr<readdy::model::actions::AddParticles>
    addParticles(const std::vector<readdy::model::Particle> &particles) const override {
        return std::make_unique<readdy::model::actions::AddParticles>(kernel, particles);
    }

    std::unique_ptr<readdy::model::actions::EulerBDIntegrator> eulerBDIntegrator(scalar timeStep) const override;

    std::unique_ptr<readdy::model::actions::MdgfrdIntegrator> mdgfrdIntegrator(scalar timeStep) const override;

    std::unique_ptr<readdy::model::actions::CalculateForces> calculateForces() const override;

    std::unique_ptr<readdy::model::actions::CreateNeighborList>
    createNeighborList(scalar interactionDistance) const override;

    std::unique_ptr<readdy::model::actions::UpdateNeighborList> updateNeighborList() const override;

    std::unique_ptr<readdy::model::actions::ClearNeighborList> clearNeighborList() const override;

    std::unique_ptr<readdy::model::actions::EvaluateCompartments> evaluateCompartments() const override;

    std::unique_ptr<readdy::model::actions::reactions::UncontrolledApproximation>
    uncontrolledApproximation(scalar timeStep) const override;

    std::unique_ptr<readdy::model::actions::reactions::Gillespie>
    gillespie(scalar timeStep) const override;

    std::unique_ptr<readdy::model::actions::reactions::DetailedBalance>
    detailedBalance(scalar timeStep) const override;

    std::unique_ptr<readdy::model::actions::top::EvaluateTopologyReactions>
    evaluateTopologyReactions(scalar timeStep) const override;

private:
    SCPUKernel *const kernel;
};
}
}
}
}
