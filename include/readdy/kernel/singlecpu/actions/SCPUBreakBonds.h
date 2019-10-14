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

namespace readdy::kernel::scpu::actions::top {

class SCPUBreakBonds : public readdy::model::actions::top::BreakBonds {
public:
    explicit SCPUBreakBonds(SCPUKernel *kernel, scalar timeStep, readdy::model::actions::top::BreakConfig config)
            : BreakBonds(timeStep, std::move(config)), kernel(kernel) {}

    void perform() override {
        std::size_t topologyIdx = 0;
        for (auto *top : kernel->stateModel().getTopologies()) {
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

                auto resultingTopologies = reaction.execute(*top, kernel);

                auto &topologies = kernel->getSCPUKernelStateModel().topologies();
                auto &model = kernel->getSCPUKernelStateModel();
                const auto &context = kernel->context();

                if (!resultingTopologies.empty()) {
                    // we had a topology fission, so we need to actually remove the current topology from the
                    // data structure
                    topologies.erase(topologies.begin() + topologyIdx);
                    //log::error("erased topology with index {}", event.topology_idx);
                    assert(top->isDeactivated());

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

                } else {
                    if (top->isNormalParticle(*kernel)) {
                        kernel->getSCPUKernelStateModel().getParticleData()->entry_at(
                                top->getParticles().front()).topology_index = -1;
                        topologies.erase(topologies.begin() + topologyIdx);
                        //log::error("erased topology with index {}", event.topology_idx);
                        assert(top->isDeactivated());
                    }
                }
            }
            ++topologyIdx;
        }
    }

private:
    SCPUKernel *kernel;

    // Note that this also evaluates forces as a by-product, i.e. adds to the force property of particles
    scalar
    evaluateEdgeEnergy(std::tuple<vertex_ref, vertex_ref> edge, const readdy::model::top::GraphTopology &t) const {
        const auto[vertex1, vertex2] = edge;

        // find bond configurations for given edge
        std::unordered_map<api::BondType, std::vector<readdy::model::top::pot::BondConfiguration>, readdy::util::hash::EnumClassHash> bondConfigs;
        {
            const auto &potentialConfiguration = kernel->context().topologyRegistry().potentialConfiguration();
            auto it = potentialConfiguration.pairPotentials.find(
                    std::tie(vertex1->particleType(), vertex2->particleType()));
            if (it != potentialConfiguration.pairPotentials.end()) {
                for (const auto &cfg : it->second) {
                    bondConfigs[cfg.type].emplace_back(vertex1->particleIndex, vertex2->particleIndex,
                                                       cfg.forceConstant, cfg.length);
                }
            } else {
                std::ostringstream ss;
                auto p1 = t.particleForVertex(vertex1);
                auto p2 = t.particleForVertex(vertex2);

                ss << "The edge " << vertex1->particleIndex << " (" << kernel->context().particleTypes().nameOf(p1.type()) << ")";
                ss << " -- " << vertex2->particleIndex << " (" << kernel->context().particleTypes().nameOf(p2.type()) << ")";
                ss << " has no bond configured!";

                throw std::invalid_argument(ss.str());
            }
        }

        // transform configurations to potential instances
        std::vector<std::unique_ptr<readdy::model::top::Topology::bonded_potential>> bondedPotentials;
        for (const auto &bond : bondConfigs) {
            switch (bond.first) {
                case api::BondType::HARMONIC: {
                    bondedPotentials.push_back(std::make_unique<readdy::model::top::TopologyActionFactory::harmonic_bond>(bond.second));
                    break;
                };
            }
        }

        // create actions, perform and accumulate energies on edge
        auto taf = kernel->getTopologyActionFactory();
        scalar totalEnergyForEdge {0.};
        for (const auto &bondedPot : bondedPotentials) {
            totalEnergyForEdge += bondedPot->createForceAndEnergyAction(taf)->perform(&t);
        }

        return totalEnergyForEdge;
    }
};

}
