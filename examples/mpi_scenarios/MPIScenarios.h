/**
 * Scenarios to measure performance, and optionally observables to ensure that the physics is correct
 * Examples to implement with constant density and adaptive volume based on number of workers:
 * - MSD -> free diffusion
 * - Stationary distribution in an external potential -> diffusion in double well along one dimension
 * - Thermodynamics of LJ suspension -> diffusion subject to pair-interactions
 * - Michaelis-Menten kinetics -> irreversible reaction-diffusion without forces in the reaction-limit (well-mixed)
 * - A+B<-->C with LJ stationary distr. -> reversible reaction-diffusion with forces, can do diffusion-influenced
 *
 * @file Scenarios.h
 * @brief << brief description >>
 * @author chrisfroe
 * @date 28.02.20
 */

#pragma once

#include <readdy/model/Kernel.h>
#include <readdy/kernel/mpi/MPIKernel.h>
#include <readdy/api/Simulation.h>
#include <readdy/api/KernelConfiguration.h>
#include <utility>
#include <readdy/kernel/mpi/model/MPIUtils.h>
#include <Scenarios.h>

namespace readdy::kernel::mpi::benchmark {

using Json = nlohmann::json;

class MPIDistributeParticles : public readdy::performance::Scenario {
public:
    MPIDistributeParticles() : Scenario(
            "MPIDistributeParticles",
            "Distribute particles") {}

    Json run() override {
        int worldSize;
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        // find out before context so we can adjust the load (number of particles) according to weak scaling setup
        std::size_t nWorkers = worldSize - 1;
        std::size_t nParticlesPerWorker = 10000;

        readdy::model::Context ctx;

        ctx.boxSize() = {static_cast<readdy::scalar>(nWorkers * 5.), 5., 5.};
        ctx.particleTypes().add("A", 1.);
        ctx.particleTypes().add("B", 1.);
        ctx.potentials().addHarmonicRepulsion("A", "A", 10., 2.4);
        Json conf = {{"MPI", {{"dx", 4.9}, {"dy", 4.9}, {"dz", 4.9}}}};
        ctx.kernelConfiguration() = conf.get<readdy::conf::Configuration>();

        readdy::kernel::mpi::MPIKernel kernel(ctx);

        assert(nWorkers == kernel.domain().nWorkerRanks());

        auto idA = kernel.context().particleTypes().idOf("A");
        const std::size_t nParticles = nParticlesPerWorker * nWorkers;
        std::vector<readdy::model::Particle> particles;
        for (std::size_t i = 0; i < nParticles; ++i) {
            auto x = readdy::model::rnd::uniform_real() * ctx.boxSize()[0] - 0.5 * ctx.boxSize()[0];
            auto y = readdy::model::rnd::uniform_real() * 5. - 2.5;
            auto z = readdy::model::rnd::uniform_real() * 5. - 2.5;
            particles.emplace_back(x, y, z, idA);
        }

        auto addParticles = kernel.actions().addParticles(particles);
        MPI_Barrier(kernel.commUsedRanks());
        {
            readdy::util::Timer t("addParticles");
            addParticles->perform();
        }

        if (kernel.domain().isMasterRank()) {
            assert(kernel.getMPIKernelStateModel().getParticleData()->size() == 0); // master data is emtpy
        } else if (kernel.domain().isWorkerRank()) {
            assert(kernel.getMPIKernelStateModel().getParticleData()->size() > 0); // worker should have gotten one particle
        } else if (kernel.domain().isIdleRank()) {
            assert(kernel.getMPIKernelStateModel().getParticleData()->size() == 0); // idle workers are idle
        } else {
            throw std::runtime_error("Must be one of those above");
        }

        const auto currentParticles = kernel.getMPIKernelStateModel().gatherParticles();
        if (kernel.domain().isMasterRank()) {
            assert(currentParticles.size() == nParticles);
        }

        Json result;
        result["context"] = ctx.describe();
        result["performance"] = Json::parse(readdy::util::Timer::perfToJsonString());
        readdy::util::Timer::clear();
        return result;
    }
};

}
