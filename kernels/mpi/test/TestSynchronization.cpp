#include <catch2/catch.hpp>
#include <readdy/kernel/mpi/MPIKernel.h>
#include <readdy/testing/Utils.h>

/**
 * << detailed description >>
 *
 * @file TestActions.cpp
 * @brief << brief description >>
 * @author chrisfroe
 * @date 04.03.20
 */

using json = nlohmann::json;

void check(const readdy::kernel::mpi::MPIKernel &kernel, std::size_t expNumberTotal, std::size_t expNumberResp) {
    const auto data = kernel.getMPIKernelStateModel().getParticleData();
    if (kernel.domain()->isWorkerRank()) {
        auto n = std::count_if(data->begin(), data->end(),
                               [](const auto &entry) { return !entry.deactivated; });
        CHECK(n == expNumberTotal);
        auto nResp = std::count_if(data->begin(), data->end(),
                                   [](const auto &entry) { return !entry.deactivated and entry.responsible; });
        CHECK(nResp == expNumberResp);
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
}

using ParticlePODSet = std::unordered_set<readdy::kernel::mpi::util::ParticlePOD, readdy::kernel::mpi::util::HashPOD>;

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
                {0.5,  0., 0., idA}, // is in halo of neighbor
                {1.5,  0., 0., idA},
                {2.5,  0., 0., idA},
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

        WHEN("States are synchronized") {
            kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
            THEN("each rank should see 4 particles, 3 for which it is responsible and 1 in halo") {
                check(kernel, 4, 3);
            }
            AND_WHEN("States are synchronized again") {
                kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
                THEN("We see the same result") {
                    check(kernel, 4, 3);
                }
            }
            AND_WHEN("Particles are gathered again") {
                auto gatheredParticles = kernel.getMPIKernelStateModel().gatherParticles();
                THEN("All 6 particles are obtained") {
                    if (kernel.domain()->isMasterRank()) {
                        CHECK(gatheredParticles.size() == 6);
                    }
                }
            }
        }
    }

    GIVEN("Three domains in one periodic dimension with 15 particles") {
        readdy::model::Context ctx;
        ctx.boxSize() = {15., 5., 5.};
        ctx.periodicBoundaryConditions() = {true, false, false};
        json conf = {{"MPI", {{"dx", 4.9}, {"dy", 4.9}, {"dz", 4.9}}}};
        ctx.kernelConfiguration() = conf.get<readdy::conf::Configuration>();
        ctx.particleTypes().add("A", 1.);
        ctx.potentials().addHarmonicRepulsion("A", "A", 1., 1.);
        readdy::kernel::mpi::MPIKernel kernel(ctx);

        auto idA = kernel.context().particleTypes().idOf("A");
        std::vector<readdy::model::Particle> particles = {
                {-7, 0., 0., idA}, // is in halo of neighbor
                {-6, 0., 0., idA},
                {-5, 0., 0., idA},
                {-4, 0., 0., idA},
                {-3, 0., 0., idA}, // is in halo of neighbor
                // domain boundary
                {-2, 0., 0., idA}, // is in halo of neighbor
                {-1, 0., 0., idA},
                {0,  0., 0., idA},
                {1,  0., 0., idA},
                {2,  0., 0., idA}, // is in halo of neighbor
                // domain boundary
                {3,  0., 0., idA}, // is in halo of neighbor
                {4,  0., 0., idA},
                {5,  0., 0., idA},
                {6,  0., 0., idA},
                {7,  0., 0., idA}, // is in halo of neighbor
        };
        kernel.getMPIKernelStateModel().distributeParticles(particles);

        WHEN("Particles are gathered") {
            auto gatheredParticles = kernel.getMPIKernelStateModel().gatherParticles();
            THEN("All 15 particles are obtained") {
                if (kernel.domain()->isMasterRank()) {
                    CHECK(gatheredParticles.size() == 15);
                }
            }
        }

        WHEN("States are synchronized") {
            kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
            THEN("Each worker sees 5 responsible particles and 2 halo particles") {
                check(kernel, 7, 5);
            }

            AND_WHEN("States are synchronized again") {
                kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
                THEN("Each worker still sees 5 responsible particles and 2 halo particles") {
                    check(kernel, 7, 5);
                }
            }

            AND_WHEN("Particles are gathered again") {
                auto gatheredParticles = kernel.getMPIKernelStateModel().gatherParticles();
                THEN("All 15 particles are obtained") {
                    if (kernel.domain()->isMasterRank()) {
                        CHECK(gatheredParticles.size() == 15);
                    }
                }
            }

            THEN("Worker 3 has the correct particles at the correct position") {
                if (kernel.domain()->rank() == 3) {
                    ParticlePODSet expected;
                    expected.emplace(particles[0]);
                    expected.emplace(particles[9]);
                    expected.emplace(particles[10]);
                    expected.emplace(particles[11]);
                    expected.emplace(particles[12]);
                    expected.emplace(particles[13]);
                    expected.emplace(particles[14]);

                    ParticlePODSet actual;
                    for (const auto &entry : *(kernel.getMPIKernelStateModel().getParticleData())) {
                        if (!entry.deactivated) {
                            actual.emplace(entry);
                        }
                    }

                    CHECK(expected == actual);
                }

            }
        }
    }
}
