/**
 * Scenarios to measure performance (profiling, scaling),
 * and optionally observables to ensure that the physics is correct.
 *
 * Scenarios to ensure correct physics could be e.g.:
 * - MSD -> free diffusion
 * - Stationary distribution in an external potential -> diffusion in double well along one dimension
 * - Thermodynamics of LJ suspension -> diffusion subject to pair-interactions
 * - Michaelis-Menten kinetics -> irreversible reaction-diffusion without forces in the reaction-limit (well-mixed)
 * - A+B<-->C with LJ stationary distr. -> reversible reaction-diffusion with forces, can do diffusion-influenced
 *
 * @file MPIScenarios.h
 * @brief Some scenarios to run to measure performance (profiling, scaling) or observables
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

namespace readdy::performance {

using Json = nlohmann::json;

class MPIDistributeParticles : public Scenario {
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

class MPIDiffusionPairPotential : public Scenario {
    WeakScalingGeometry _mode;
    // parameter determining the volume and thus the final number of particles
    // describes the ratio of box length and interaction distance
    readdy::scalar _edgeLengthOverInteractionDistance;
public:
    explicit MPIDiffusionPairPotential(WeakScalingGeometry mode, readdy::scalar edgeLengthOverInteractionDistance = 5.) : Scenario(
            "MPIDiffusionPairPotential",
            "Diffusion of particles with constant density, scale box according to mode"), _mode(mode), _edgeLengthOverInteractionDistance(edgeLengthOverInteractionDistance) {}

    Json run() override {
        int rank, worldSize;
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::size_t nLoad = worldSize; // in weak scaling the load and processors must scale linearly
        std::size_t nProcessors = worldSize;
        std::size_t nWorkers = worldSize - 1;
        std::size_t nParticles;
        std::size_t nParticlesPerLoad;
        scalar nParticlesPerWorker;
        std::array<readdy::scalar, 3> box{};
        scalar halo = 2.;
        // radius that is used to calculate the volume occupation, if the particles were hard-spheres
        scalar assumedRadius = halo / 2.;
        if (_mode == stick) {
            // will lead to 60% volume occupation when one particle has an associated radius of halo/2.
            nParticlesPerLoad = 0.6 * std::pow(_edgeLengthOverInteractionDistance, 3)
                                / (4./3. * readdy::util::numeric::pi<scalar>() * std::pow(assumedRadius/halo,3));
            nParticles = nLoad * nParticlesPerLoad;
            nParticlesPerWorker = static_cast<scalar>(nParticles) / static_cast<scalar>(nWorkers);
            scalar boxLength = _edgeLengthOverInteractionDistance * halo;
            box = {nLoad * boxLength, boxLength, boxLength};
            readdy::log::info("rank={}, stick mode, nParticlesPerLoad {}, nParticles {}", rank, nParticlesPerLoad, nParticles);
        } else if (_mode == cube) {
            // layman's check if nWorkers is a cubic number
            auto remainder = fmod(std::cbrt(nWorkers), 1.0f);
            readdy::log::critical("rank={}, nWorkers={}, remainder={}", rank, nWorkers, remainder);
            if (std::abs(remainder) > 0.0001) {
                throw std::invalid_argument("The supplied number of workers will not entirely fill a cubic geometry");
            }
            nParticlesPerLoad = 0.6 * std::pow(_edgeLengthOverInteractionDistance, 3)
                                / (4./3. * readdy::util::numeric::pi<scalar>() * std::pow(assumedRadius/halo,3));
            nParticles = nLoad * nParticlesPerLoad;
            nParticlesPerWorker = static_cast<scalar>(nParticles) / static_cast<scalar>(nWorkers);
            scalar boxLength = _edgeLengthOverInteractionDistance * halo;
            scalar factor = std::cbrt(static_cast<scalar>(nLoad));
            box = {factor * boxLength, factor * boxLength, factor * boxLength};
            readdy::log::info("rank={}, cube mode, nParticlesPerLoad {}, nParticles {}", rank, nParticlesPerLoad, nParticles);
        } else {
            throw std::invalid_argument(fmt::format("Unknown scaling mode {}", _mode));
        }

        readdy::model::Context ctx;
        ctx.boxSize() = box;
        ctx.particleTypes().add("A", 1.);
        ctx.potentials().addHarmonicRepulsion("A", "A", 10., halo);

        readdy::kernel::mpi::MPIKernel kernel(ctx);

        assert(nLoad == kernel.domain().worldSize());

        auto idA = kernel.context().particleTypes().idOf("A");
        std::vector<readdy::model::Particle> particles;
        for (std::size_t i = 0; i < nParticles; ++i) {
            auto x = readdy::model::rnd::uniform_real() * ctx.boxSize()[0] - 0.5 * ctx.boxSize()[0];
            auto y = readdy::model::rnd::uniform_real() * ctx.boxSize()[1] - 0.5 * ctx.boxSize()[1];
            auto z = readdy::model::rnd::uniform_real() * ctx.boxSize()[2] - 0.5 * ctx.boxSize()[2];
            particles.emplace_back(x, y, z, idA);
        }

        auto addParticles = kernel.actions().addParticles(particles);
        MPI_Barrier(MPI_COMM_WORLD);
        {
            readdy::util::Timer t("addParticles");
            addParticles->perform();
        }

        readdy::scalar timeStep = 0.01;
        std::size_t nSteps = 1000;
        auto integrator = kernel.actions().eulerBDIntegrator(timeStep);
        auto forces = kernel.actions().calculateForces();
        auto neighborList = kernel.actions().updateNeighborList();

        neighborList->perform();
        forces->perform();
        MPI_Barrier(MPI_COMM_WORLD);
        for (size_t t = 1; t < nSteps + 1; t++) {
            readdy::util::Timer tStep("complete timestep");
            {
                readdy::util::Timer t1("integrator");
                integrator->perform();
            }
            {
                readdy::util::Timer t2("neighborList");
                neighborList->perform();
            }
            {
                readdy::util::Timer t3("forces");
                forces->perform();
            }
        }

        Json result;
        result["context"] = ctx.describe();
        result["domain"] = kernel.domain().describe();
        result["performance"] = Json::parse(readdy::util::Timer::perfToJsonString());
        result["nLoad"] = nLoad;
        result["nProcessors"] = nProcessors;
        result["nParticles"] = nParticles;
        result["nParticlesPerProcessor"] = static_cast<scalar>(nParticles) / static_cast<scalar>(nProcessors);
        result["nParticlesPerLoad"] = nParticlesPerLoad;
        result["nParticlesPerWorkers"] = nParticlesPerWorker;
        result["kernelName"] = "MPI";
        result["mode"] = fmt::format("{}", _mode);
        result["edgeLengthOverInteractionDistance"] = _edgeLengthOverInteractionDistance;
        readdy::util::Timer::clear();
        return result;
    }
};

}
