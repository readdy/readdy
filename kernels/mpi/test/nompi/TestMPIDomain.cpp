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
 * @file TestMPIDomain.cpp
 * @brief Test proper construction of domains and their neighborhood
 * @author chrisfroe
 * @date 17.06.19
 */

#include <catch2/catch.hpp>
#include <readdy/kernel/mpi/model/MPIDomain.h>

using NeighborType = readdy::kernel::mpi::model::MPIDomain::NeighborType;

TEST_CASE("Test domain decomposition", "[mpi]") {
    readdy::model::Context context;
    context.particleTypes().add("A", 1.0);
    context.potentials().addHarmonicRepulsion("A", "A", 1.0, 2.3); // cutoff 2.3
    
    SECTION("Test wrapIntoThisHalo") {
        // todo
    }

    SECTION("1D chain of domains") {
        context.boxSize() = {{10., 1., 1.}};
        context.periodicBoundaryConditions() = {{true, false, false}};
        // user given values should be at least twice the cutoff
        std::array<readdy::scalar, 3> userMinDomainWidths{{4.6, 4.6, 4.6}};
        // expecting two workers, boxsize will be split in half, and one master rank -> 3
        int worldSize = 3;

        for (int rank = 0; rank < worldSize; ++rank) {
            readdy::kernel::mpi::model::MPIDomain domain(rank, worldSize, userMinDomainWidths, context);
            REQUIRE(domain.worldSize == 3);
            REQUIRE(domain.nDomains() == std::array<std::size_t, 3>({{2, 1, 1}}));
            REQUIRE(domain.domainIndex()(0, 0, 0) == 0);
            REQUIRE(domain.rank == rank);
            REQUIRE(domain.haloThickness == context.calculateMaxCutoff());
            if (rank != 0) {
                readdy::Vec3 five{{4.99 - 5., 0., 0.}};
                readdy::Vec3 twoandahalf{{2.5 - 5., 0., 0.}};
                readdy::Vec3 six{{5.01 - 5., 0., 0.}};
                readdy::Vec3 sevenandahalf{{7.5 - 5., 0., 0.}};

                if (domain.isInDomainCore(five) or domain.isInDomainCore(twoandahalf)) {
                    // if either of them is in domain, both must be in domain
                    REQUIRE(domain.isInDomainCore(five));
                    REQUIRE(domain.isInDomainCore(twoandahalf));

                    // six is not in this domain but very close
                    REQUIRE(domain.isInDomainHalo(six));
                    // sevenandahalf however is too far in the next domain, and thus not in this ones' halo, nor core
                    REQUIRE_FALSE(domain.isInDomainHalo(sevenandahalf));
                    REQUIRE_FALSE(domain.isInDomainCoreOrHalo(sevenandahalf));
                    REQUIRE_FALSE(domain.isInDomainCore(sevenandahalf));
                }

                if (domain.isInDomainCore(six) or domain.isInDomainCore(sevenandahalf)) {
                    REQUIRE(domain.isInDomainCore(six));
                    REQUIRE(domain.isInDomainCore(sevenandahalf));

                    REQUIRE(domain.isInDomainHalo(five));
                    REQUIRE_FALSE(domain.isInDomainHalo(twoandahalf));
                    REQUIRE_FALSE(domain.isInDomainCoreOrHalo(twoandahalf));
                    REQUIRE_FALSE(domain.isInDomainCore(twoandahalf));
                }
            }
            // check neighborhood

            if (rank != 0) {
                auto otherRank = rank == 1 ? 2 : 1;
                // neighborIndex (2,1,1) is the neighbor in the west (x+dx, y, z)
                REQUIRE(domain.neighborRanks()[domain.neighborIndex(1 + 1, 1, 1)] == otherRank);
                // (0, 1, 1) is the neighbor in the east
                REQUIRE(domain.neighborRanks()[domain.neighborIndex(1 - 1, 1, 1)] == otherRank);
                // "neighbor" (1,1,1) is this domain/rank
                REQUIRE(domain.neighborRanks()[domain.neighborIndex(1, 1, 1)] == rank);

                REQUIRE(domain.neighborTypes()[domain.neighborIndex(1 + 1, 1, 1)] == NeighborType::regular);
                REQUIRE(domain.neighborTypes()[domain.neighborIndex(1 - 1, 1, 1)] == NeighborType::regular);
                REQUIRE(domain.neighborTypes()[domain.neighborIndex(1, 1, 1)] == NeighborType::self);

                // all other neighbors do not exist -> ranks -1, types nan
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        for (int k = 0; k < 3; ++k) {
                            // skip the ones we checked above
                            if ((not(i == 1 and j == 1 and k == 1))
                                and (not(i == 0 and j == 1 and k == 1))
                                and (not(i == 2 and j == 1 and k == 1))) {
                                REQUIRE(domain.neighborRanks()[domain.neighborIndex(i, j, k)] == -1);
                                REQUIRE(domain.neighborTypes()[domain.neighborIndex(i, j, k)] == NeighborType::nan);
                            }
                        }
                    }
                }
            }
        }

        SECTION("Rong user given domain widths") {
            std::array<readdy::scalar, 3> incompatibleDomainWidths{{2.2, 4.6, 4.6}};
            auto ctor = [&]() {
                readdy::kernel::mpi::model::MPIDomain domain(0, 4 + 1, incompatibleDomainWidths, context);
            };
            REQUIRE_THROWS(ctor());
        }
    }

    SECTION("Check one position in (4 by 2 by 1) domains, periodic in xy") {
        context.boxSize() = {{20., 10., 1.}};
        context.periodicBoundaryConditions() = {{true, true, false}};
        std::array<readdy::scalar, 3> userMinDomainWidths{{4.6, 4.6, 4.6}};
        int worldSize = 1 + (4 * 2 * 1);
        std::size_t posCount {0};
        readdy::Vec3 pos {{-8., -4., 0.49}};
        int south, north, west, east, northwest, southwest;
        for (int rank = 0; rank < worldSize; ++rank) {
            readdy::kernel::mpi::model::MPIDomain domain(rank, worldSize, userMinDomainWidths, context);
            REQUIRE(domain.worldSize == 9);
            REQUIRE(domain.nDomains() == std::array<std::size_t, 3>({{4, 2, 1}}));
            REQUIRE(domain.domainIndex()(0, 0, 0) == 0);
            REQUIRE(domain.rank == rank);
            REQUIRE(domain.haloThickness == context.calculateMaxCutoff());
            if (rank != 0) {
                if (domain.isInDomainCore(pos)) {
                    posCount++;

                    west =  domain.neighborRanks()[domain.neighborIndex(0,1,1)];
                    east =  domain.neighborRanks()[domain.neighborIndex(2,1,1)];
                    south = domain.neighborRanks()[domain.neighborIndex(1,0,1)];
                    north = domain.neighborRanks()[domain.neighborIndex(1,2,1)];
                    northwest =  domain.neighborRanks()[domain.neighborIndex(0,2,1)];
                    southwest =  domain.neighborRanks()[domain.neighborIndex(0,0,1)];

                    REQUIRE(northwest == southwest); // due to periodicity
                    REQUIRE(north == south); // due to periodicity
                }
            }
        }
        REQUIRE(posCount == 1); // can only be in one domain core

        // todo specify if global or local position

        readdy::kernel::mpi::model::MPIDomain southDomain(south, worldSize, userMinDomainWidths, context);
        REQUIRE(southDomain.isInDomainHalo(pos));
        readdy::kernel::mpi::model::MPIDomain northDomain(north, worldSize, userMinDomainWidths, context);
        REQUIRE(northDomain.isInDomainHalo(pos));
        readdy::kernel::mpi::model::MPIDomain westDomain(west, worldSize, userMinDomainWidths, context);
        REQUIRE(westDomain.isInDomainHalo(pos));
        readdy::kernel::mpi::model::MPIDomain eastDomain(east, worldSize, userMinDomainWidths, context);
        REQUIRE_FALSE(eastDomain.isInDomainHalo(pos));
        readdy::kernel::mpi::model::MPIDomain nwDomain(northwest, worldSize, userMinDomainWidths, context);
        REQUIRE(nwDomain.isInDomainHalo(pos));
    }

    SECTION("Periodic in a direction (here y) which has only one domain") {
        context.boxSize() = {{10., 5., 1.}};
        context.periodicBoundaryConditions() = {{true, true, false}};
        // user given values should be at least twice the cutoff
        std::array<readdy::scalar, 3> userMinDomainWidths{{4.6, 4.6, 4.6}};
        int worldSize = 3;
        for (int rank = 0; rank < worldSize; ++rank) {
            readdy::kernel::mpi::model::MPIDomain domain(rank, worldSize, userMinDomainWidths, context);
            REQUIRE(domain.worldSize == 3);
            REQUIRE(domain.nDomains() == std::array<std::size_t, 3>({{2, 1, 1}}));
            REQUIRE(domain.domainIndex()(0, 0, 0) == 0);
            REQUIRE(domain.rank == rank);
            REQUIRE(domain.haloThickness == context.calculateMaxCutoff());
            if (rank != 0) {
                readdy::Vec3 five{{4.99 - 5., -0.1, 0.}};
                readdy::Vec3 twoandahalf{{2.5 - 5., -0.1, 0.}};
                readdy::Vec3 six{{5.01 - 5., -0.1, 0.}};
                readdy::Vec3 sevenandahalf{{7.5 - 5., -0.1, 0.}};

                if (domain.isInDomainCore(five) or domain.isInDomainCore(twoandahalf)) {
                    // if either of them is in domain, both must be in domain
                    REQUIRE(domain.isInDomainCore(five));
                    REQUIRE(domain.isInDomainCore(twoandahalf));

                    // six is not in this domain but very close
                    REQUIRE(domain.isInDomainHalo(six));
                    // sevenandahalf however is too far in the next domain, and thus not in this ones' halo, nor core
                    REQUIRE_FALSE(domain.isInDomainHalo(sevenandahalf));
                    REQUIRE_FALSE(domain.isInDomainCoreOrHalo(sevenandahalf));
                    REQUIRE_FALSE(domain.isInDomainCore(sevenandahalf));
                }

                if (domain.isInDomainCore(six) or domain.isInDomainCore(sevenandahalf)) {
                    REQUIRE(domain.isInDomainCore(six));
                    REQUIRE(domain.isInDomainCore(sevenandahalf));

                    REQUIRE(domain.isInDomainHalo(five));
                    REQUIRE_FALSE(domain.isInDomainHalo(twoandahalf));
                    REQUIRE_FALSE(domain.isInDomainCoreOrHalo(twoandahalf));
                    REQUIRE_FALSE(domain.isInDomainCore(twoandahalf));
                }
            }
            // check neighborhood
            if (rank != 0) {
                auto otherRank = rank == 1 ? 2 : 1;

                auto west = otherRank; // two boxes and periodic in WE direction
                auto east = otherRank;
                auto south = rank; // only one box and periodic in SN direction
                auto north = rank;

                REQUIRE(domain.neighborRanks()[domain.neighborIndex(1 - 1, 1, 1)] == west);
                REQUIRE(domain.neighborRanks()[domain.neighborIndex(1 + 1, 1, 1)] == east);
                REQUIRE(domain.neighborRanks()[domain.neighborIndex(1, 1 - 1, 1)] == south);
                REQUIRE(domain.neighborRanks()[domain.neighborIndex(1, 1 + 1, 1)] == north);

                // up and down are not periodic and very thin -> only one domain on this axis
                REQUIRE(domain.neighborRanks()[domain.neighborIndex(1, 1, 1 - 1)] == -1);
                REQUIRE(domain.neighborRanks()[domain.neighborIndex(1, 1, 1 + 1)] == -1);


                REQUIRE(domain.neighborRanks()[domain.neighborIndex(1, 1, 1)] == rank);

                REQUIRE(domain.neighborTypes()[domain.neighborIndex(1 + 1, 1, 1)] == NeighborType::regular);
                REQUIRE(domain.neighborTypes()[domain.neighborIndex(1 - 1, 1, 1)] == NeighborType::regular);
                REQUIRE(domain.neighborTypes()[domain.neighborIndex(1, 1 - 1, 1)] == NeighborType::self);
                REQUIRE(domain.neighborTypes()[domain.neighborIndex(1, 1 + 1, 1)] == NeighborType::self);
                REQUIRE(domain.neighborTypes()[domain.neighborIndex(1, 1, 1 - 1)] == NeighborType::nan);
                REQUIRE(domain.neighborTypes()[domain.neighborIndex(1, 1, 1 + 1)] == NeighborType::nan);
                REQUIRE(domain.neighborTypes()[domain.neighborIndex(1, 1, 1)] == NeighborType::self);

                // requesting a position outside the box makes no sense, except either along west or east
                auto center = domain.origin() + 0.5 * domain.extent();
                readdy::Vec3 dEW{{3., 0., 0.}};
                readdy::Vec3 dNS{{0., 3., 0.}};
                readdy::Vec3 dUD{{0., 0., 3.}};
                REQUIRE_THROWS(domain.rankOfPosition(center - dUD));
                REQUIRE_THROWS(domain.rankOfPosition(center + dUD));
                REQUIRE_THROWS(domain.rankOfPosition(center - dNS));
                REQUIRE_THROWS(domain.rankOfPosition(center + dNS));
                if (center[0] < 0.) {
                    REQUIRE_THROWS(domain.rankOfPosition(center - dEW));
                    REQUIRE(domain.rankOfPosition(center + dEW) == otherRank);
                } else {
                    REQUIRE_THROWS(domain.rankOfPosition(center + dEW));
                    REQUIRE(domain.rankOfPosition(center - dEW) == otherRank);
                }
            }
        }
    }
}
