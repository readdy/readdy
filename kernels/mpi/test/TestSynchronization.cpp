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
        WHEN("particles are distributed and gathered again we ")
        kernel.getMPIKernelStateModel().distributeParticles(particles);
        // todo assert that particles are at given locations
        kernel.getMPIKernelStateModel().gatherParticles();

        kernel.getMPIKernelStateModel().synchronizeWithNeighbors();
        // todo assert that domain1 sees 4 particles (3 own + 1 in halo), same for domain2
    }
}

TEST_CASE("Test send and receive operations", "[mpi]") {

}