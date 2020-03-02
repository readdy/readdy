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
 * A general complication with running MPI tests is that the number of processes is set from the outside via mpirun,
 * so the tests should be called with a certain minimum amount of workers. Having more workers is OK, some will
 * just idle. Having less will cause an error because all domains must have a worker.
 *
 * @file TestSetupMPI.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 28.05.19
 */

#include <catch2/catch.hpp>

#include <readdy/model/Kernel.h>
#include <readdy/kernel/mpi/MPIKernel.h>
#include <readdy/api/Simulation.h>
#include <readdy/api/KernelConfiguration.h>

using json = nlohmann::json;

TEST_CASE("Kernel construction and delayed initialization", "[mpi]") {
    WHEN("Kernel is constructed without context") {
        readdy::kernel::mpi::MPIKernel kernel;
        THEN("The kernel is not initialized") {
            CHECK_FALSE(kernel.isInitialized());
        }
    }
}

TEST_CASE("Test mpi kernel observe particle number", "[mpi]") {
    MPI_Barrier(MPI_COMM_WORLD);
    if (false) {
        readdy::model::Context ctx;
        readdy::Simulation simulation("MPI", ctx);

        REQUIRE(simulation.selectedKernelType() == "MPI");

        SECTION("In and out types and positions") {
            ctx.boxSize() = {10., 10., 10.};
            ctx.particleTypes().add("A", 1.);
            ctx.particleTypes().add("B", 1.);
            //ctx.reactions().add("fusili: A +(1.) A -> B", 0.1);
            ctx.potentials().addHarmonicRepulsion("A", "A", 10., 2.3);
            json conf = {{"MPI", {{"dx", 4.9}, {"dy", 4.9}, {"dz", 4.9}}}};
            ctx.kernelConfiguration() = conf.get<readdy::conf::Configuration>();
            const std::size_t nParticles = 100;
            for (std::size_t i = 0; i < nParticles; ++i) {
                auto x = readdy::model::rnd::uniform_real() * 10. - 5.;
                auto y = readdy::model::rnd::uniform_real() * 10. - 5.;
                auto z = readdy::model::rnd::uniform_real() * 10. - 5.;
                simulation.addParticle("A", x, y, z);
            }
            if (false) {
                const auto idA = ctx.particleTypes().idOf("A");
                auto check = [&nParticles, &idA](readdy::model::observables::Particles::result_type result) {
                    const auto &types = std::get<0>(result);
                    const auto &ids = std::get<1>(result);
                    const auto &positions = std::get<2>(result);
                    CHECK(std::count(types.begin(), types.end(), idA));
                };
                auto obsHandle = simulation.registerObservable(simulation.observe().particles(1), check);
                simulation.run(100, 0.01);
            }
        }
    }
}

TEST_CASE("Test kernel configuration from context", "[mpi]") {
    MPI_Barrier(MPI_COMM_WORLD);
    readdy::model::Context ctx;
    ctx.boxSize() = {10., 10., 10.};
    ctx.particleTypes().add("A", 0.1);
    ctx.particleTypes().add("B", 0.1);
    ctx.potentials().addHarmonicRepulsion("B", "B", 1.0, 2.);
    json conf = {{"MPI", {{"dx", 4.9}, {"dy", 4.9}, {"dz", 4.9}}}};
    ctx.kernelConfiguration() = conf.get<readdy::conf::Configuration>();

    CHECK(ctx.kernelConfiguration().mpi.dx == Approx(4.9));
    CHECK(ctx.kernelConfiguration().mpi.dy == Approx(4.9));
    CHECK(ctx.kernelConfiguration().mpi.dz == Approx(4.9));
    readdy::kernel::mpi::MPIKernel kernel(ctx); // this also initializes domains

    CHECK(kernel.context().kernelConfiguration().mpi.dx == Approx(4.9));
    CHECK(kernel.context().kernelConfiguration().mpi.dy == Approx(4.9));
    CHECK(kernel.context().kernelConfiguration().mpi.dz == Approx(4.9));
}

TEST_CASE("Test distribute particles and gather them again", "[mpi]") {
    MPI_Barrier(MPI_COMM_WORLD);
    readdy::model::Context ctx;

    ctx.boxSize() = {10., 10., 10.};
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    ctx.potentials().addHarmonicRepulsion("A", "A", 10., 2.3);
    json conf = {{"MPI", {{"dx", 4.9}, {"dy", 4.9}, {"dz", 4.9}}}};
    ctx.kernelConfiguration() = conf.get<readdy::conf::Configuration>();

    readdy::kernel::mpi::MPIKernel kernel(ctx); // this also initializes domains

    auto idA = kernel.context().particleTypes().idOf("A");
    const std::size_t nParticles = kernel.domain()->nDomains(); // one particle per domain
    for (const auto &rank : kernel.domain()->workerRanks()) {
        auto [origin, extent] = kernel.domain()->coreOfDomain(rank);
        auto center = origin + 0.5 * extent;
        readdy::model::Particle p(center, idA);
        kernel.getMPIKernelStateModel().addParticle(p);
    }

    if (kernel.domain()->isMasterRank()) {
        CHECK(kernel.getMPIKernelStateModel().getParticleData()->size() == 0); // master data is emtpy
    } else if (kernel.domain()->isWorkerRank()) {
        CHECK(kernel.getMPIKernelStateModel().getParticleData()->size() == 1); // worker should have gotten one particle
    } else if (kernel.domain()->isIdleRank()) {
        CHECK(kernel.getMPIKernelStateModel().getParticleData()->size() == 0); // idle workers are idle
    } else {
        throw std::runtime_error("Must be one of those above");
    }

    const auto currentParticles = kernel.stateModel().getParticles();

    if (kernel.domain()->isMasterRank()) {
        CHECK(currentParticles.size() == nParticles);
        for (const auto &rank : kernel.domain()->workerRanks()) {
            readdy::Vec3 origin, extent;
            std::tie(origin, extent) = kernel.domain()->coreOfDomain(rank);
            auto center = origin + 0.5 * extent;
            readdy::model::Particle p(center, idA);
            auto found = std::find_if(currentParticles.begin(), currentParticles.end(),
                                      [&p](const readdy::model::Particle &p2) -> bool {
                                          return p.type() == p2.type() and (p.pos() - p2.pos()).normSquared() < 0.01;
                                      });
            CHECK(found != currentParticles.end());
        }
    }
}
