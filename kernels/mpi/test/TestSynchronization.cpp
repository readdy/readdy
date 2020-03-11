#include <catch2/catch.hpp>
#include <readdy/kernel/mpi/MPIKernel.h>
#include <readdy/testing/Utils.h>
#include <readdy/common/boundary_condition_operations.h>

/**
 * << detailed description >>
 *
 * @file TestActions.cpp
 * @brief << brief description >>
 * @author chrisfroe
 * @date 04.03.20
 */

namespace rkm = readdy::kernel::mpi;
namespace rkmu = readdy::kernel::mpi::util;

struct HashPODPair {
    std::size_t operator()(const std::pair<rkmu::ParticlePOD,rkmu::ParticlePOD> &podPair) const {
        auto [p1, p2] = podPair;
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

using Json = nlohmann::json;
using ParticlePODSet = std::unordered_set<rkmu::ParticlePOD, rkmu::HashPOD>;
using ParticlePODPairSet = std::unordered_set<std::pair<rkmu::ParticlePOD,rkmu::ParticlePOD>, HashPODPair>;

void check(const readdy::kernel::mpi::MPIKernel &kernel, std::size_t expNumberTotal, std::size_t expNumberResp,
           const ParticlePODSet &expectedPODs = {}, const ParticlePODPairSet &expectedPODPairs = {}) {
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

void synchronizeAndCheck(readdy::kernel::mpi::MPIKernel &kernel, std::size_t expNumberTotal, std::size_t expNumberResp,
                         const ParticlePODSet &expectedPODs, const ParticlePODPairSet &expectedPODPairs,
                         const std::vector<readdy::model::Particle> &allParticles) {
    WHEN("Particles are gathered again") {
        auto gatheredParticles = kernel.getMPIKernelStateModel().gatherParticles();
        THEN("All 6 particles are obtained") {
            if (kernel.domain()->isMasterRank()) {
                CHECK(gatheredParticles.size() == allParticles.size());
            }
        }
    }

    WHEN("States are synchronized") {
        kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
        THEN("each rank should see nPerWorkerTotal particles, nPerWorkerResp for which it is responsible") {
            check(kernel, expNumberTotal, expNumberResp);
        }
        AND_WHEN("States are synchronized again") {
            kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
            THEN("We see the same result") {
                check(kernel, expNumberTotal, expNumberResp);
            }
        }
        AND_WHEN("Particles are gathered again") {
            auto gatheredParticles = kernel.getMPIKernelStateModel().gatherParticles();
            THEN("All particles are obtained") {
                if (kernel.domain()->isMasterRank()) {
                    CHECK(gatheredParticles.size() == allParticles.size());
                }
            }
        }
        THEN("Worker has the correct particles at the correct position, and the correct pairs") {
            if (kernel.domain()->isWorkerRank()) {
                ParticlePODSet actual;
                for (const auto &entry : *(kernel.getMPIKernelStateModel().getParticleData())) {
                    if (!entry.deactivated) {
                        actual.emplace(entry);
                    }
                }
                CHECK(expectedPODs == actual);

                // todo pairs
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
        for (std::size_t i=0; i<particles.size(); ++i) {
            const auto& pi = particles.at(i);

            // for single particles, the worker needs to see all that are in domain core or halo
            if (kernel.domain()->isInDomainCoreOrHalo(pi.pos())) {
                expectedPODs.emplace(pi);
            }

            // for pairs we only consider core-core and core-halo pairs but not halo-halo,
            // because the worker will not evaluate forces/reactions for those
            if (kernel.domain()->isInDomainCore(pi.pos())) {
                for (std::size_t j=i+1; j<particles.size(); ++j) {
                    const auto& pj = particles.at(j);
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
        std::size_t nPerWorkerTotal = 4; // todo extract this from the set, have another set for responsible particles
        std::size_t nPerWorkerResp = 3;

        auto [expectedPODs, expectedPODPairs] = expectedParticlesAndPairs(kernel, particles);

        kernel.getMPIKernelStateModel().distributeParticles(particles);

        synchronizeAndCheck(kernel, nPerWorkerTotal, nPerWorkerResp, expectedPODs, expectedPODPairs, particles);
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

        std::size_t nPerWorkerTotal = 7;
        std::size_t nPerWorkerResp = 5;

        auto [expectedPODs, expectedPODPairs] = expectedParticlesAndPairs(kernel, particles);

        kernel.getMPIKernelStateModel().distributeParticles(particles);

        synchronizeAndCheck(kernel, nPerWorkerTotal, nPerWorkerResp, expectedPODs, expectedPODPairs, particles);
    }

    // todo more cases with multiple domains in 2D and 3D, exclude certain cases should something be raised?
    // todo how about 1 domain in a coordinate with periodic boundary?
}
