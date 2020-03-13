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
 * << detailed description >>
 *
 * @file TestNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 1/18/18
 */

#include <catch2/catch.hpp>

#include <readdy/model/Context.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>
#include <readdy/kernel/singlecpu/model/topologies/SCPUTopologyActionFactory.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/common/boundary_condition_operations.h>

TEST_CASE("Test the singlecpu neighbor list", "[scpu]") {
    using namespace readdy;
    using IndexPair = std::tuple<std::size_t, std::size_t>;

    auto boxSize = GENERATE(
            std::array<readdy::scalar, 3>{{15, 15, 15}},
            std::array<readdy::scalar, 3>{{10, 10, 10}},
            std::array<readdy::scalar, 3>{{10, 5, 5}},
            std::array<readdy::scalar, 3>{{5, 5, 10}},
            std::array<readdy::scalar, 3>{{15, 5, 10}},
            std::array<readdy::scalar, 3>{{5, 10, 5}},
            std::array<readdy::scalar, 3>{{5, 5, 5}});

    INFO(fmt::format("Testing for box size ({}, {}, {})", boxSize[0], boxSize[1], boxSize[2]));

    kernel::scpu::SCPUKernel kernel {};
    auto& context = kernel.context();
    context.particleTypes().add("Test", 1.);
    auto id = context.particleTypes().idOf("Test");
    scalar cutoff = 4;
    context.reactions().addFusion("Fusion", id, id, id, .001, cutoff);
    context.boxSize() = boxSize;

    bool periodic = true;
    context.periodicBoundaryConditions()[0] = periodic;
    context.periodicBoundaryConditions()[1] = periodic;
    context.periodicBoundaryConditions()[2] = periodic;

    kernel::scpu::model::top::SCPUTopologyActionFactory taf (nullptr);

#ifdef READDY_DEBUG
    auto n_steps = 3U;
    auto n_particles = 1000;
#else
    auto n_steps = 5U;
    auto n_particles = 3000;
#endif

    for(auto i = 0; i < n_particles; ++i) {
        model::Particle particle(model::rnd::uniform_real<scalar>(-.5*context.boxSize()[0], .5*context.boxSize()[0]),
                                 model::rnd::uniform_real<scalar>(-.5*context.boxSize()[1], .5*context.boxSize()[1]),
                                 model::rnd::uniform_real<scalar>(-.5*context.boxSize()[2], .5*context.boxSize()[2]), id);
        auto& stateModel = kernel.getSCPUKernelStateModel();
        stateModel.addParticle(particle);
    }

    kernel.initialize();

    kernel.stateModel().initializeNeighborList(context.calculateMaxCutoff());

    auto integrator = kernel.actions().eulerBDIntegrator(.1);
    auto reactionHandler = kernel.actions().uncontrolledApproximation(.1);

    const auto &data = *kernel.getSCPUKernelStateModel().getParticleData();

    // collect all pairs of particles that are closer than cutoff, these should (uniquely) be in the NL
    std::unordered_set<IndexPair,
            util::ForwardBackwardTupleHasher<IndexPair>,
            util::ForwardBackwardTupleEquality<IndexPair>> pairs;
    for(auto t = 0U; t < n_steps; ++t) {
        integrator->perform();

        kernel.stateModel().updateNeighborList();
        {
            pairs.clear();
            std::size_t ix1 = 0;
            for(const auto &e1 : data) {
                std::size_t ix2 = 0;
                for(const auto &e2 : data) {
                    if(ix1 != ix2 && !e1.deactivated && !e2.deactivated &&
                       bcs::dist(e1.pos, e2.pos, context.boxSize().data(), context.periodicBoundaryConditions().data()) < cutoff) {
                        pairs.insert(std::make_tuple(ix1, ix2));
                    }
                    ++ix2;
                }
                ++ix1;
            }
        }
        auto pairsCopy = pairs;

        const auto &neighborList = *kernel.getSCPUKernelStateModel().getNeighborList();
        for (auto cell = 0U; cell < neighborList.nCells(); ++cell) {
            for(auto it = neighborList.particlesBegin(cell); it != neighborList.particlesEnd(cell); ++it) {
                const auto &entry = data.entry_at(*it);
                // A deactivated entry should not end up in the NL
                REQUIRE_FALSE(entry.deactivated);
                neighborList.forEachNeighbor(it, cell, [&](const std::size_t neighborIndex) {
                    const auto &neighbor = data.entry_at(neighborIndex);
                    // A deactivated entry should not end up in the NL
                    REQUIRE_FALSE(neighbor.deactivated);
                    // we got a pair
                    if(bcs::dist(entry.pos, neighbor.pos, context.boxSize().data(), context.periodicBoundaryConditions().data()) < cutoff) {
                        auto findIt = pairs.find(std::make_tuple(*it, neighborIndex));
                        REQUIRE(findIt != pairs.end());
                        if(findIt != pairs.end()) {
                            pairs.erase(findIt);
                        }
                    }
                });
            }
        }
        // all pairs should end up in the NL
        REQUIRE(pairs.empty());

        reactionHandler->perform();

    }
}
