/********************************************************************
 * Copyright © 2020 Noe Group, Freie Universität Berlin (GER)       *
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

#include <catch2/catch.hpp>
#include <readdy/kernel/mpi/MPIKernel.h>

/**
 * @file TestDiffusion.cpp
 * @brief Simulation loop with particles diffusing subject to forces
 * @author chrisfroe
 * @date 13.03.20
 */

namespace rkmu = readdy::kernel::mpi::util;
namespace rkm = readdy::kernel::mpi;
namespace rnd = readdy::model::rnd;

TEST_CASE("Test diffusion conservation of particles when diffusing", "[mpi]") {
    GIVEN("System of A and B particles subject to soft repulsion") {
        readdy::model::Context ctx;
        ctx.boxSize() = {10., 10., 10.};
        ctx.periodicBoundaryConditions() = {true, true, true};
        ctx.particleTypes().add("A", 1.0);
        ctx.particleTypes().add("B", 1.0);
        ctx.potentials().addHarmonicRepulsion("A", "B", 1.0, 1.0);

        rkm::MPIKernel kernel(ctx);

        auto idA = kernel.context().particleTypes().idOf("A");
        auto idB = kernel.context().particleTypes().idOf("B");
        const auto &box = kernel.context().boxSize();
        std::vector<readdy::model::Particle> particles;
        std::size_t na{100}, nb{200};

        for (std::size_t i = 0; i < na; ++i) {
            readdy::Vec3 pos{rnd::uniform_real() * box[0] - 0.5 * box[0],
                             rnd::uniform_real() * box[1] - 0.5 * box[1],
                             rnd::uniform_real() * box[2] - 0.5 * box[2]};
            particles.emplace_back(pos, idA);
        }
        for (std::size_t i = 0; i < nb; ++i) {
            readdy::Vec3 pos{rnd::uniform_real() * box[0] - 0.5 * box[0],
                             rnd::uniform_real() * box[1] - 0.5 * box[1],
                             rnd::uniform_real() * box[2] - 0.5 * box[2]};
            particles.emplace_back(pos, idB);
        }

        WHEN("Particles diffuse") {
            readdy::scalar timeStep = 0.01;
            auto integrator = kernel.actions().eulerBDIntegrator(timeStep);
            auto forces = kernel.actions().calculateForces();
            auto neighborList = kernel.actions().updateNeighborList();
            // adding particles as part of an Action to separate collective operations from the user
            auto addParticles = kernel.actions().addParticles(particles);
            addParticles->perform();

            std::size_t nSteps = 1000;

            neighborList->perform();
            forces->perform();
            // todo observables -> evaluate results and gather(0) on workers, gather(0) on master
            //kernel.evaluateObservables(0); todo
            for (size_t t = 1; t < nSteps + 1; t++) {
                integrator->perform();
                neighborList->perform();
                forces->perform();
                //kernel.evaluateObservables(t); todo
            }
            THEN("Number of particles is conserved") {
                auto ps = kernel.getMPIKernelStateModel().gatherParticles();
                if (kernel.domain().isMasterRank()) {
                    auto numberA = std::count_if(ps.begin(), ps.end(), [&](const readdy::model::Particle& p) {return p.type() == idA;});
                    auto numberB = std::count_if(ps.begin(), ps.end(), [&](const readdy::model::Particle& p) {return p.type() == idB;});
                    CHECK(numberA == na);
                    CHECK(numberB == nb);
                }
            }
        }
    }
}
