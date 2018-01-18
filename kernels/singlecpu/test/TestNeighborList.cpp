/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * @file TestNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 1/18/18
 */

#include <gtest/gtest.h>
#include <readdy/model/Context.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>
#include <readdy/kernel/singlecpu/model/topologies/SCPUTopologyActionFactory.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace {

TEST(NeighborList, Diffuse) {
    using namespace readdy;

    using IndexPair = std::tuple<std::size_t, std::size_t>;

    kernel::scpu::SCPUKernel kernel;

    auto& context = kernel.context();
    context.particle_types().add("Test", 1.);
    auto id = context.particle_types().idOf("Test");
    scalar cutoff = 1;
    context.reactions().addFusion("Fusion", id, id, id, .001, cutoff);
    context.boxSize()[0] = 10;
    context.boxSize()[1] = 10;
    context.boxSize()[2] = 10;
    bool periodic = true;
    context.periodicBoundaryConditions()[0] = periodic;
    context.periodicBoundaryConditions()[1] = periodic;
    context.periodicBoundaryConditions()[2] = periodic;

    kernel::scpu::model::top::SCPUTopologyActionFactory taf (nullptr);

#ifdef READDY_DEBUG
    auto n_steps = 10U;
    auto n_particles = 1000;
#else
    auto n_steps = 20U;
    auto n_particles = 3000;
#endif

    for(auto i = 0; i < n_particles; ++i) {
        model::Particle particle(model::rnd::uniform_real<scalar>(-5, 5),
                                 model::rnd::uniform_real<scalar>(-5, 5),
                                 model::rnd::uniform_real<scalar>(-5, 5), id);
        kernel.stateModel().addParticle(particle);
    }

    kernel.initialize();

    kernel.stateModel().initializeNeighborList(0);

    auto integrator = kernel.actions().eulerBDIntegrator(.1);
    auto reactionHandler = kernel.actions().uncontrolledApproximation(.1);

    const auto &d2 = kernel.context().distSquaredFun();
    const auto &data = *kernel.getSCPUKernelStateModel().getParticleData();

    for(auto t = 0U; t < n_steps; ++t) {
        integrator->perform();
        
        kernel.stateModel().updateNeighborList();
        // collect all pairs of particles that are closer than cutoff, these should (uniquely) be in the NL
        std::unordered_set<IndexPair,
                           util::ForwardBackwardTupleHasher<IndexPair>,
                           util::ForwardBackwardTupleEquality<IndexPair>> pairs;
        {
            std::size_t ix1 = 0;
            for(const auto &e1 : data) {
                std::size_t ix2 = 0;
                for(const auto &e2 : data) {
                    if(ix1 != ix2 && !e1.deactivated && !e2.deactivated && std::sqrt(d2(e1.pos, e2.pos)) < cutoff) {
                        pairs.insert(std::make_tuple(ix1, ix2));
                    }
                    ++ix2;
                }
                ++ix1;
            }
        }

        const auto &neighborList = *kernel.getSCPUKernelStateModel().getNeighborList();
        for (auto cell = 0U; cell < neighborList.nCells(); ++cell) {
            for(auto it = neighborList.particlesBegin(cell); it != neighborList.particlesEnd(cell); ++it) {
                const auto &entry = data.entry_at(*it);
                EXPECT_FALSE(entry.deactivated) << "A deactivated entry should not end up in the NL";
                neighborList.forEachNeighbor(it, cell, [&](const std::size_t neighborIndex) {
                    const auto &neighbor = data.entry_at(neighborIndex);
                    EXPECT_FALSE(neighbor.deactivated) << "A deactivated entry should not end up in the NL";
                    // we got a pair
                    if(std::sqrt(d2(entry.pos, neighbor.pos)) < cutoff) {
                        auto findIt = pairs.find(std::make_tuple(*it, neighborIndex));
                        EXPECT_NE(findIt, pairs.end()) << "A particle pairing was more than once in the NL";
                        if(findIt != pairs.end()) {
                            pairs.erase(findIt);
                        }
                    }
                });
            }
        }
        EXPECT_EQ(pairs.size(), 0) << "Some pairs were not contained in the NL";

        reactionHandler->perform();

        /*std::cout << "n_particles: " << std::count_if(data.begin(), data.end(), [](const auto &entry) {
            return !entry.deactivated;
        }) << std::endl;*/
    }
}

}