#include <catch2/catch.hpp>
#include <readdy/kernel/mpi/MPIKernel.h>

/**
 * << detailed description >>
 *
 * @file TestActions.cpp
 * @brief << brief description >>
 * @author chrisfroe
 * @date 04.03.20
 */

using json = nlohmann::json;

TEST_CASE("Synchronization of neighbors", "[mpi]") {
    GIVEN("Two domains with halo-thickness=1 and 6 particles, two of which are in the other domain's halo") {
        readdy::model::Context ctx;
        ctx.boxSize() = {10., 5., 5.};
        ctx.periodicBoundaryConditions() = {false, false, false};
        json conf = {{"MPI", {{"dx", 4.9}, {"dy", 4.9}, {"dz", 4.9}}}};
        ctx.kernelConfiguration() = conf.get<readdy::conf::Configuration>();
        ctx.particleTypes().add("A", 1.);
        ctx.potentials().addHarmonicRepulsion("A", "A", 1., 1.);
        readdy::kernel::mpi::MPIKernel kernel(ctx);

        auto idA = kernel.context().particleTypes().idOf("A");
        std::vector<readdy::model::Particle> particles = {
                {-2.5, 0., 0., idA},
                {-1.5, 0., 0., idA},
                {-0.5, 0., 0., idA}, // is in halo of neighbor
                {0.5, 0., 0., idA}, // is in halo of neighbor
                {1.5, 0., 0., idA},
                {2.5, 0., 0., idA},
        };
        kernel.getMPIKernelStateModel().distributeParticles(particles);

        WHEN("Particles are gathered again") {
            auto gatheredParticles = kernel.getMPIKernelStateModel().gatherParticles();
            THEN("All 6 particles are obtained") {
                if (kernel.domain()->isMasterRank()) {
                    CHECK(gatheredParticles.size() == 6);
                }
            }
        }

        const auto check = [&]() {
            const auto data = kernel.getMPIKernelStateModel().getParticleData();
            if (kernel.domain()->isWorkerRank()) {
                auto n = std::count_if(data->begin(), data->end(),
                                       [](const auto &entry) { return !entry.deactivated; });
                CHECK(n == 4);
                auto nResp = std::count_if(data->begin(), data->end(), [](const auto &entry){ return !entry.deactivated and entry.responsible; });
                CHECK(nResp == 3);
                for (const auto &entry : *data) {
                    if (!entry.deactivated) {
                        if (entry.responsible) {
                            CHECK(kernel.domain()->isInDomainCore(entry.position()));
                        } else {
                            CHECK(kernel.domain()->isInDomainHalo(entry.position()));
                        }
                    }
                }
            } else if (kernel.domain()->isMasterRank()) {
                // master has no active entries
                CHECK(data->size() == data->n_deactivated());
            }
        };

        WHEN("States are synchronized") {
            kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
            THEN("each rank should see 4 particles, 3 for which it is responsible and 1 in halo") {
                check();
            }
            AND_WHEN("States are synchronized again") {
                kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
                THEN("We see the same result") {
                    check();
                }
            }
        }
    }
}

TEST_CASE("Test send and receive operations", "[mpi]") {

}