
/**
 * << detailed description >>
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

namespace readdy::kernel::mpi::performance {

using Json = nlohmann::json;

// There is no proper testing here, just some sanity requirements
void require(bool statement, std::string description = "") {
    if (!statement) {
        throw std::logic_error("Statement is false, " + description);
    }
}

class Scenario {
    std::string _description{};
    std::string _name{};
public:
    explicit Scenario(std::string name, std::string descr) : _name(std::move(name)), _description(std::move(descr)) {

    }

    const std::string &description() {
        return _description;
    };

    const std::string &name() {
        return _name;
    }

    virtual Json run() = 0;
};

class DistributeParticles : public Scenario {
public:
    DistributeParticles() : Scenario(
            "DistributeParticles",
            "Distribute particles and gather them again") {}

    Json run() override {
        int worldSize;
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        std::size_t nWorkers = worldSize - 1;
        std::size_t nParticlesPerWorker = 10000;

        readdy::model::Context ctx;

        ctx.boxSize() = {static_cast<readdy::scalar>(nWorkers * 5.), 5., 5.};
        ctx.particleTypes().add("A", 1.);
        ctx.particleTypes().add("B", 1.);
        ctx.potentials().addHarmonicRepulsion("A", "A", 10., 2.4);
        Json conf = {{"MPI", {{"dx", 4.9}, {"dy", 4.9}, {"dz", 4.9}}}};
        ctx.kernelConfiguration() = conf.get<readdy::conf::Configuration>();

        readdy::kernel::mpi::MPIKernel kernel(ctx); // this also initializes domains

        require(kernel.domain()->nWorkerRanks() == nWorkers);
        require(kernel.domain()->worldSize() == worldSize);

        auto idA = kernel.context().particleTypes().idOf("A");
        const std::size_t nParticles = nParticlesPerWorker * nWorkers;
        std::vector<readdy::model::Particle> particles;
        for (std::size_t i = 0; i < nParticles; ++i) {
            auto x = readdy::model::rnd::uniform_real() * ctx.boxSize()[0] - 0.5 * ctx.boxSize()[0];
            auto y = readdy::model::rnd::uniform_real() * 5. - 2.5;
            auto z = readdy::model::rnd::uniform_real() * 5. - 2.5;
            particles.emplace_back(x, y, z, idA);
        }

        MPI_Barrier(kernel.commUsedRanks());
        kernel.stateModel().addParticles(particles);

        if (kernel.domain()->isMasterRank()) {
            require(kernel.getMPIKernelStateModel().getParticleData()->size() == 0); // master data is emtpy
        } else if (kernel.domain()->isWorkerRank()) {
            require(kernel.getMPIKernelStateModel().getParticleData()->size() >
                    0); // worker should have gotten one particle
        } else if (kernel.domain()->isIdleRank()) {
            require(kernel.getMPIKernelStateModel().getParticleData()->size() == 0); // idle workers are idle
        } else {
            throw std::runtime_error("Must be one of those above");
        }

        const auto currentParticles = kernel.stateModel().getParticles();
        if (kernel.domain()->isMasterRank()) {
            require(currentParticles.size() == nParticles);
        }

        Json json = readdy::kernel::mpi::util::Timer::perfToJson();
        readdy::kernel::mpi::util::Timer::clear();
        return json;
    }
};

}
