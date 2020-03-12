#include <catch2/catch.hpp>
#include <readdy/kernel/mpi/MPIKernel.h>
#include <readdy/common/boundary_condition_operations.h>

/**
 * Test correct distribution of particles from master to workers, and synchronization of particles between workers
 *
 * @file TestSynchronization.cpp
 * @brief Test synchronization of particles from master to workers and in between workers
 * @author chrisfroe
 * @date 04.03.20
 */

namespace rkm = readdy::kernel::mpi;
namespace rkmu = readdy::kernel::mpi::util;
namespace rnd = readdy::model::rnd;
using Json = nlohmann::json;

struct HashPODPair {
    std::size_t operator()(const std::pair<rkmu::ParticlePOD, rkmu::ParticlePOD> &podPair) const {
        auto[p1, p2] = podPair;
        std::size_t s1 = rkmu::HashPOD{}(p1);
        std::size_t s2 = rkmu::HashPOD{}(p2);
        std::size_t seed{0};
        // permutation invariant hash
        if (s1 > s2) {
            readdy::util::hash::combine(seed, s1);
            readdy::util::hash::combine(seed, s2);
        } else {
            readdy::util::hash::combine(seed, s2);
            readdy::util::hash::combine(seed, s1);
        }
        return seed;
    }
};

struct ComparePODPair {
    bool operator()(
            const std::pair<rkmu::ParticlePOD, rkmu::ParticlePOD> &podPair1,
            const std::pair<rkmu::ParticlePOD, rkmu::ParticlePOD> &podPair2) {
        std::size_t s1 = HashPODPair{}(podPair1);
        std::size_t s2 = HashPODPair{}(podPair2);
        return s1 < s2;
    }
};

using ParticlePODSet = std::unordered_set<rkmu::ParticlePOD, rkmu::HashPOD>;
using ParticlePODPairSet = std::unordered_set<std::pair<rkmu::ParticlePOD, rkmu::ParticlePOD>, HashPODPair>;

void check(readdy::kernel::mpi::MPIKernel &kernel,
           const ParticlePODSet &expectedPODs, const ParticlePODPairSet &expectedPODPairs) {
    const auto data = kernel.getMPIKernelStateModel().getParticleData();
    if (kernel.domain()->isWorkerRank()) {
        // cannot use sections inside diverging parts of the program, thus commented out
        // THEN("The number of particles is correct and responsible particles are in the domain core")
        auto n = std::count_if(data->begin(), data->end(),
                               [](const auto &entry) { return !entry.deactivated; });
        CHECK(n == expectedPODs.size());
        for (const auto &entry : *data) {
            if (!entry.deactivated) {
                if (entry.responsible) {
                    CHECK(kernel.domain()->isInDomainCore(entry.position()));
                } else {
                    CHECK(kernel.domain()->isInDomainHalo(entry.position()));
                }
            }
        }
        // THEN("Worker has the correct particles at the correct position")
        {
            ParticlePODSet actual;
            for (const auto &entry : *(kernel.getMPIKernelStateModel().getParticleData())) {
                if (!entry.deactivated) {
                    actual.emplace(entry);
                }
            }
            CHECK(expectedPODs == actual);
        }
        // AND_WHEN("Neighborlist is filled")
        {
            readdy::scalar willBeIgnored = 1.;
            kernel.getMPIKernelStateModel().initializeNeighborList(willBeIgnored);

            //THEN("the correct pairs are in the neighborlist")
            {
                ParticlePODPairSet actual;

                auto emplacePair = [&actual](rkm::MPIEntry e1, rkm::MPIEntry e2) {
                    actual.emplace(e1, e2);
                };

                const auto &nl = kernel.getMPIKernelStateModel().getNeighborList();
                nl->forAllPairs(emplacePair);

                std::vector<std::pair<rkmu::ParticlePOD, rkmu::ParticlePOD>> actualVec(actual.begin(), actual.end());
                std::vector<std::pair<rkmu::ParticlePOD, rkmu::ParticlePOD>> expectedVec(expectedPODPairs.begin(),
                                                                                         expectedPODPairs.end());

                std::sort(actualVec.begin(), actualVec.end(), ComparePODPair{});
                std::sort(expectedVec.begin(), expectedVec.end(), ComparePODPair{});
                // expected should be a subset of actual
                CHECK(std::includes(actualVec.begin(), actualVec.end(), expectedVec.begin(), expectedVec.end(), ComparePODPair{}));
            }
        }
    } else if (kernel.domain()->isMasterRank()) {
        // master has no active entries
        CHECK(data->size() == data->n_deactivated());
    }
}

void synchronizeAndCheck(readdy::kernel::mpi::MPIKernel &kernel,
                         const ParticlePODSet &expectedPODs, const ParticlePODPairSet &expectedPODPairs,
                         const std::vector<readdy::model::Particle> &allParticles) {
    WHEN("Particles are gathered") {
        auto gatheredParticles = kernel.getMPIKernelStateModel().gatherParticles();
        THEN("All particles are obtained") {
            if (kernel.domain()->isMasterRank()) {
                CHECK(gatheredParticles.size() == allParticles.size());
            }
        }
    }
    WHEN("States are synchronized") {
        kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
        THEN("The state is correctly set up") {
            check(kernel, expectedPODs, expectedPODPairs);
            AND_WHEN("States are synchronized again") {
                kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
                THEN("We see the same result") {
                    check(kernel, expectedPODs, expectedPODPairs);
                    AND_WHEN("Particles are gathered again") {
                        auto gatheredParticles = kernel.getMPIKernelStateModel().gatherParticles();
                        THEN("All particles are obtained") {
                            if (kernel.domain()->isMasterRank()) {
                                CHECK(gatheredParticles.size() == allParticles.size());
                            }
                        }
                    }
                }
            }
        }
    }
}

std::pair<ParticlePODSet, ParticlePODPairSet> expectedParticlesAndPairs(
        const readdy::kernel::mpi::MPIKernel &kernel,
        const std::vector<readdy::model::Particle> &particles) {
    if (kernel.domain()->isWorkerRank()) {
        ParticlePODSet expectedPODs;
        ParticlePODPairSet minimalExpectedPODPairs; // the workers have to have these at least
        for (std::size_t i = 0; i < particles.size(); ++i) {
            const auto &pi = particles.at(i);

            // for single particles, the worker needs to see all that are in domain core or halo
            if (kernel.domain()->isInDomainCoreOrHalo(pi.pos())) {
                expectedPODs.emplace(pi);
            }

            // for pairs we only consider core-core and core-halo pairs but not halo-halo,
            // because the worker will not evaluate forces/reactions for those
            if (kernel.domain()->isInDomainCore(pi.pos())) {
                for (std::size_t j = i + 1; j < particles.size(); ++j) {
                    const auto &pj = particles.at(j);
                    if (kernel.domain()->isInDomainCoreOrHalo(pj.pos())) {
                        // assume that the only interaction has a radius of 1
                        auto distance = readdy::bcs::dist(
                                pi.pos(), pj.pos(), kernel.context().boxSize(),
                                kernel.context().periodicBoundaryConditions());
                        if (distance <= 1.0) {
                            minimalExpectedPODPairs.emplace(pi, pj);
                        }
                    }
                }
            }
        }
        return {expectedPODs, minimalExpectedPODPairs};
    } else {
        return {};
    }
}

void setupContext(readdy::model::Context &ctx) {
    Json conf = {{"MPI", {{"dx", 4.9}, {"dy", 4.9}, {"dz", 4.9}}}};
    ctx.kernelConfiguration() = conf.get<readdy::conf::Configuration>();
    ctx.particleTypes().add("A", 1.);
    ctx.potentials().addHarmonicRepulsion("A", "A", 1., 1.);
}

TEST_CASE("Synchronization of neighbors", "[mpi]") {
    GIVEN("Two domains with halo-thickness=1 and 6 particles, two of which are in the other domain's halo") {
        readdy::model::Context ctx;

        ctx.boxSize() = {10., 5., 5.};
        ctx.periodicBoundaryConditions() = {false, false, false};

        setupContext(ctx);
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

        auto[expectedPODs, expectedPODPairs] = expectedParticlesAndPairs(kernel, particles);
        kernel.getMPIKernelStateModel().distributeParticles(particles);
        synchronizeAndCheck(kernel, expectedPODs, expectedPODPairs, particles);
    }

    GIVEN("Three domains in one periodic dimension with 15 particles") {
        readdy::model::Context ctx;

        ctx.boxSize() = {15., 5., 5.};
        ctx.periodicBoundaryConditions() = {true, false, false};

        setupContext(ctx);
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
        auto[expectedPODs, expectedPODPairs] = expectedParticlesAndPairs(kernel, particles);
        kernel.getMPIKernelStateModel().distributeParticles(particles);
        synchronizeAndCheck(kernel, expectedPODs, expectedPODPairs, particles);
    }

    GIVEN("Three domains with all periodic dimensions") {
        readdy::model::Context ctx;

        ctx.boxSize() = {15., 5., 5.};
        ctx.periodicBoundaryConditions() = {true, true, true};

        setupContext(ctx);
        readdy::kernel::mpi::MPIKernel kernel(ctx);
        auto idA = kernel.context().particleTypes().idOf("A");

        std::vector<readdy::model::Particle> particles = {
                {-7, 0., 0., idA},
                {-6, 0., 0., idA},
                {-5, 0., 0., idA},
                {-4, 0., 0., idA},
                {-3, 0., 0., idA},
                {-2, 0., 0., idA},
                {-1, 0., 0., idA},
                {0,  0., 0., idA},
                {0,  0., 2.4, idA},
                {1,  0., 0., idA},
                {2,  0., 0., idA},
                {3,  0., 0., idA},
                {4,  0., 0., idA},
                {5,  0., 0., idA},
                {6,  0., 0., idA},
                {7,  0., 0., idA},
        };
        auto[expectedPODs, expectedPODPairs] = expectedParticlesAndPairs(kernel, particles);
        kernel.getMPIKernelStateModel().distributeParticles(particles);
        synchronizeAndCheck(kernel, expectedPODs, expectedPODPairs, particles);
    }
}

TEST_CASE("Two dimensional synchronization", "[mpi]") {
    readdy::model::Context ctx;

    ctx.boxSize() = {15.01, 15.01, 5.};
    ctx.periodicBoundaryConditions() = {false, false, false};

    setupContext(ctx);
    readdy::kernel::mpi::MPIKernel kernel(ctx);
    auto idA = kernel.context().particleTypes().idOf("A");

    std::vector<readdy::model::Particle> particles = {
            // one line of particles along x
            {-2, 0, 0, idA},
            {-1, 0, 0, idA},
            {0,  0, 0, idA},
            {1,  0, 0, idA},
            {2,  0, 0, idA},
            // one line of particles along y
            {0, -2, 0, idA},
            {0, -1, 0, idA},
            {0,  0.01, 0, idA},
            {0,  1, 0, idA},
            {0,  2, 0, idA},
    };

    auto[expectedPODs, expectedPODPairs] = expectedParticlesAndPairs(kernel, particles);
    kernel.getMPIKernelStateModel().distributeParticles(particles);
    synchronizeAndCheck(kernel, expectedPODs, expectedPODPairs, particles);
}

void randomPosititionsTest(std::array<readdy::scalar, 3> boxSize, std::array<bool, 3> pbc) {
    readdy::model::Context ctx;

    ctx.boxSize() = boxSize;
    ctx.periodicBoundaryConditions() = pbc;

    setupContext(ctx);
    readdy::kernel::mpi::MPIKernel kernel(ctx);
    auto idA = kernel.context().particleTypes().idOf("A");

    const auto &box = kernel.context().boxSize();
    std::vector<readdy::model::Particle> particles;
    std::vector<rkmu::ParticlePOD> buffer;
    std::size_t nParticles = 1000;
    if (kernel.domain()->isMasterRank()) {
        for (std::size_t i=0; i<nParticles; ++i) {
            readdy::Vec3 pos{rnd::uniform_real() * box[0] - 0.5 * box[0],
                             rnd::uniform_real() * box[1] - 0.5 * box[1],
                             rnd::uniform_real() * box[2] - 0.5 * box[2]};
            particles.emplace_back(pos, idA);
            buffer.emplace_back(pos, idA);
        }
        MPI_Bcast(buffer.data(), nParticles * sizeof(rkmu::ParticlePOD), MPI_BYTE, 0, kernel.commUsedRanks());
    } else if (kernel.domain()->isWorkerRank()) {
        buffer.resize(nParticles);
        MPI_Bcast(buffer.data(), nParticles * sizeof(rkmu::ParticlePOD), MPI_BYTE, 0, kernel.commUsedRanks());
        for (const auto& pod : buffer) {
            particles.emplace_back(pod.position, pod.typeId);
        }
    } else {
        // nuffin
    }

    auto[expectedPODs, expectedPODPairs] = expectedParticlesAndPairs(kernel, particles);
    kernel.getMPIKernelStateModel().distributeParticles(particles);
    synchronizeAndCheck(kernel, expectedPODs, expectedPODPairs, particles);
}

TEST_CASE("Synchronization with random positions (have to be broadcasted) for varying box dimensions and periodicity") {
    std::vector<std::array<readdy::scalar, 3>> bs = {
            {5., 5., 5.}, // 1 worker
            {10., 5., 5.}, // 2 workers
            {10., 10., 5.}, // 4 workers
            {10., 10., 10.}, // 8 workers
            {15., 5., 10.}, // 6 workers
            {15., 5., 5.}, // 3 workers
            {15., 15., 5.}, // 9 workers
            //{15., 15., 10.}, // 18 workers
            //{15., 15., 15.}, // 27 workers
            {40., 5., 5.}, // 8 workers
    };
    std::vector<std::array<bool, 3>> pbcs = {
            {true, true, true},
            {true, true, false},
            {true, false, false},
            {true, false, true},
            {false, false, true},
            {false, true, true},
            {false, false, false},
            {false, true, false},
    };

    for (auto b : bs) {
        for (auto pbc : pbcs) {
            randomPosititionsTest(b, pbc);
        }
    }
}
