/********************************************************************
 * Copyright © 2020 Computational Molecular Biology Group,          *
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
 * @file TestObservables.cpp
 * @brief Test observables for the MPI kernel
 * @author chrisfroe
 * @date 22.04.20
 */

#include <catch2/catch.hpp>
#include <readdy/kernel/mpi/MPIKernel.h>
#include <readdy/api/Simulation.h>

using Json = nlohmann::json;

TEST_CASE("Test particles observable", "[mpi]") {
    readdy::model::Context ctx;

    /// In and out types and positions
    ctx.boxSize() = {10., 10., 10.};
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    //ctx.reactions().add("fusili: A +(1.) A -> B", 0.1);
    ctx.potentials().addHarmonicRepulsion("A", "A", 10., 2.3);
    Json conf = {{"MPI", {{"dx", 4.9}, {"dy", 4.9}, {"dz", 4.9}}}};
    ctx.kernelConfiguration() = conf.get<readdy::conf::Configuration>();

    readdy::plugin::KernelProvider::kernel_ptr kernelPtr(readdy::kernel::mpi::MPIKernel::create(ctx));
    readdy::Simulation simulation(std::move(kernelPtr));

    REQUIRE(simulation.selectedKernelType() == "MPI");

    const std::size_t nParticles = 10;
    for (std::size_t i = 0; i < nParticles; ++i) {
        auto x = readdy::model::rnd::uniform_real() * 10. - 5.;
        auto y = readdy::model::rnd::uniform_real() * 10. - 5.;
        auto z = readdy::model::rnd::uniform_real() * 10. - 5.;
        simulation.addParticle("A", x, y, z);
    }
    const auto idA = ctx.particleTypes().idOf("A");
    auto check = [&nParticles, &idA](const readdy::model::observables::Particles::result_type &result) {
        const auto &types = std::get<0>(result);
        const auto &ids = std::get<1>(result);
        const auto &positions = std::get<2>(result);
        CHECK(std::count(types.begin(), types.end(), idA) == 10);
    };
    simulation.registerObservable(simulation.observe().particles(1, check));
    simulation.run(3, 0.01);
}

// todo more tests!
