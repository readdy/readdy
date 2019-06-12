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
 * Test the general workflow of MPI kernel, i.e. set up a simulation with a couple of particles and run it.
 *
 * @file TestSetupMPI.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 28.05.19
 */

#include <catch2/catch.hpp>

#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/kernel/mpi/MPIKernel.h>
#include <readdy/api/Simulation.h>
#include <readdy/api/KernelConfiguration.h>

using json = nlohmann::json;

TEST_CASE("Test mpi kernel running in parallel", "[mpi]") {
    readdy::Simulation simulation("MPI");
    auto &ctx = simulation.context();

    REQUIRE(simulation.selectedKernelType() == "MPI");

    SECTION("In and out types and positions") {
        ctx.boxSize() = {{50.,50.,50.}};
        ctx.particleTypes().add("A", 1.);
        ctx.particleTypes().add("B", 1.);
        //ctx.reactions().add("fusili: A +(1.) A -> B", 0.1);
        ctx.potentials().addHarmonicRepulsion("A", "A", 10., 1.);
        // todo kernel configuration into context
        ctx.kernelConfiguration() = json{{"MPI", {"dx", 2.}, {"dy", 2.}, {"dz", 2.}}};
        const std::size_t nParticles = 10000;
        for (std::size_t i = 0; i < nParticles; ++i) {
            auto x = readdy::model::rnd::uniform_real()*50. - 25.;
            auto y = readdy::model::rnd::uniform_real()*50. - 25.;
            auto z = readdy::model::rnd::uniform_real()*50. - 25.;
            simulation.addParticle("A", x, y, z);
        }
        const auto idA = ctx.particleTypes().idOf("A");
        auto check = [&nParticles, &idA](readdy::model::observables::Particles::result_type result) {
            const auto &types = std::get<0>(result);
            const auto &ids = std::get<1>(result);
            const auto &positions = std::get<2>(result);
            REQUIRE(std::count(types.begin(), types.end(), idA));
        };
        auto obsHandle = simulation.registerObservable(simulation.observe().particles(1), check);
        simulation.run(100, 0.01);

    }

}