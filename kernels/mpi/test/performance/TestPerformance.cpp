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
 * « detailed description »
 *
 * @file TestPerformance.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 28.02.20
 */

#include <catch2/catch.hpp>

#include <readdy/model/Kernel.h>
#include <readdy/kernel/mpi/MPIKernel.h>
#include <readdy/api/Simulation.h>
#include <readdy/api/KernelConfiguration.h>

using json = nlohmann::json;

TEST_CASE("Test weak scaling distribute particles and gather them again ", "[!hide][mpi]") {
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    std::size_t nWorkers = worldSize - 1;
    std::size_t nParticlesPerWorker = 10000;

    readdy::model::Context ctx;

    ctx.boxSize() = {static_cast<readdy::scalar>(nWorkers * 5.), 5., 5.};
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    ctx.potentials().addHarmonicRepulsion("A", "A", 10., 2.4);
    json conf = {{"MPI", {{"dx", 4.9}, {"dy", 4.9}, {"dz", 4.9}}}};
    ctx.kernelConfiguration() = conf.get<readdy::conf::Configuration>();

    readdy::kernel::mpi::MPIKernel kernel(ctx); // this also initializes domains

    CHECK(kernel.domain()->nWorkerRanks() == nWorkers);
    CHECK(kernel.domain()->worldSize() == worldSize);

    auto idA = kernel.context().particleTypes().idOf("A");
    const std::size_t nParticles = nParticlesPerWorker * nWorkers;
    std::vector<readdy::model::Particle> particles;
    for (std::size_t i = 0; i < nParticles; ++i) {
        auto x = readdy::model::rnd::uniform_real() * ctx.boxSize()[0] - 0.5 * ctx.boxSize()[0];
        auto y = readdy::model::rnd::uniform_real() * 5. - 2.5;
        auto z = readdy::model::rnd::uniform_real() * 5. - 2.5;
        particles.emplace_back(x,y,z, idA);
    }
    kernel.stateModel().addParticles(particles);

    if (kernel.domain()->isMasterRank()) {
        CHECK(kernel.getMPIKernelStateModel().getParticleData()->size() == 0); // master data is emtpy
    } else if (kernel.domain()->isWorkerRank()) {
        CHECK(kernel.getMPIKernelStateModel().getParticleData()->size() > 0); // worker should have gotten one particle
    } else if (kernel.domain()->isIdleRank()) {
        CHECK(kernel.getMPIKernelStateModel().getParticleData()->size() == 0); // idle workers are idle
    } else {
        throw std::runtime_error("Must be one of those above");
    }

    const auto currentParticles = kernel.stateModel().getParticles();
    if (kernel.domain()->isMasterRank()) {
        CHECK(currentParticles.size() == nParticles);
    }
}
