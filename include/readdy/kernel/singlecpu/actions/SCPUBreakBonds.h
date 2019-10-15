/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
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
 * @file SCPUBreakBonds.h
 * @brief Single CPU kernel implementation of the action BreakBonds
 * @author chrisfroe
 * @date 04.10.19
 */

#pragma once

#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/model/actions/Utils.h>

namespace readdy::kernel::scpu::actions::top {

class SCPUBreakBonds : public readdy::model::actions::top::BreakBonds {
public:
    explicit SCPUBreakBonds(SCPUKernel *kernel, scalar timeStep, readdy::model::actions::top::BreakConfig config)
            : BreakBonds(timeStep, std::move(config)), kernel(kernel) {}

    void perform() override {
        auto &topologies = kernel->getSCPUKernelStateModel().topologies();
        std::vector<SCPUStateModel::topology> resultingTopologies;
        std::size_t topologyIdx = 0;
        for (auto &top : topologies) {
            if (!top->isDeactivated()) {
                auto reactionFunction = [&](
                        readdy::model::top::GraphTopology &t) -> readdy::model::top::reactions::Recipe {
                    readdy::model::top::reactions::Recipe recipe(t);
                    for (const auto &edge : t.graph().edges()) {
                        auto energy = evaluateEdgeEnergy(edge, t);
                        const auto &v1Type = std::get<0>(edge)->particleType();
                        const auto &v2Type = std::get<1>(edge)->particleType();
                        const auto &typePair = std::make_tuple(v1Type, v2Type);
                        const auto thresholdEnergyIt = thresholdEnergies().find(typePair);
                        if (thresholdEnergyIt != thresholdEnergies().end()) {
                            if (energy > thresholdEnergyIt->second) {
                                const auto &rate = breakRates().at(typePair);
                                if (readdy::model::rnd::uniform_real() < 1 - std::exp(-rate * _timeStep)) {
                                    recipe.removeEdge(edge);
                                }
                            }
                        }
                    }
                    return std::move(recipe);
                };
                scalar rateDoesntMatter{1.};
                readdy::model::top::reactions::StructuralTopologyReaction reaction("__internal_break_bonds",
                                                                                   reactionFunction, rateDoesntMatter);
                readdy::model::actions::top::executeStructuralReaction(topologies, resultingTopologies, top, reaction,
                                                                       topologyIdx,
                                                                       *(kernel->getSCPUKernelStateModel().getParticleData()),
                                                                       kernel);

            }
            ++topologyIdx;
        }

        auto &model = kernel->getSCPUKernelStateModel();
        const auto &context = kernel->context();
        for (auto &&newTopology : resultingTopologies) {
            if (!newTopology.isNormalParticle(*kernel)) {
                // we have a new topology here, update data accordingly.
                newTopology.updateReactionRates(
                        context.topologyRegistry().structuralReactionsOf(newTopology.type()));
                newTopology.configure();
                model.insert_topology(std::move(newTopology));
            } else {
                // if we have a single particle that is not of flavor topology, remove from topology structure!
                model.getParticleData()->entry_at(newTopology.getParticles().front()).topology_index = -1;
            }
        }
    }

private:
    SCPUKernel *kernel;
};

}
