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

    SECTION("1D chain of domains") {
        context.boxSize() = {10., 1., 1.};
        context.periodicBoundaryConditions() = {true, false, false};
        // user given values should be at least twice the cutoff
        context.kernelConfiguration().mpi.dx = 4.6;
        context.kernelConfiguration().mpi.dy = 4.6;
        context.kernelConfiguration().mpi.dz = 4.6;

        // expecting two workers, boxsize will be split in half, and one master rank -> 3
        int worldSize = 3;

        for (int rank = 0; rank < worldSize; ++rank) {
            context.kernelConfiguration().mpi.rank = rank;
            context.kernelConfiguration().mpi.worldSize = worldSize;
            readdy::kernel::mpi::model::MPIDomain domain(context);
            CHECK(domain.worldSize() == 3);
            CHECK(domain.nDomainsPerAxis() == std::array<std::size_t, 3>({2, 1, 1}));
            CHECK(domain.domainIndex()(0, 0, 0) == 0);
            CHECK(domain.rank() == rank);
            CHECK(domain.haloThickness() == context.calculateMaxCutoff());
            CHECK(domain.nUsedRanks() == 3);
            CHECK(domain.nWorkerRanks() == 2);
            CHECK(domain.worldSize() == domain.nUsedRanks() + domain.nIdleRanks());
            CHECK(domain.worldSize() == domain.nWorkerRanks() + 1 + domain.nIdleRanks());
            for (const auto otherRank : domain.workerRanks()) {
                readdy::Vec3 origin, extent;
                std::tie(origin, extent) = domain.coreOfDomain(otherRank);
                auto center = origin + 0.5 * extent;
                CHECK(domain.rankOfPosition(center) == otherRank);
            }
            if (rank != 0) {
                readdy::Vec3 five{4.99 - 5., 0., 0.};
                readdy::Vec3 twoandahalf{2.5 - 5., 0., 0.};
                readdy::Vec3 six{5.01 - 5., 0., 0.};
                readdy::Vec3 sevenandahalf{7.5 - 5., 0., 0.};

                if (domain.isInDomainCore(five) or domain.isInDomainCore(twoandahalf)) {
                    // if either of them is in domain, both must be in domain
                    CHECK(domain.isInDomainCore(five));
                    CHECK(domain.isInDomainCore(twoandahalf));

                    // six is not in this domain but very close
                    CHECK(domain.isInDomainHalo(six));
                    // sevenandahalf however is too far in the next domain, and thus not in this ones' halo, nor core
                    CHECK_FALSE(domain.isInDomainHalo(sevenandahalf));
                    CHECK_FALSE(domain.isInDomainCoreOrHalo(sevenandahalf));
                    CHECK_FALSE(domain.isInDomainCore(sevenandahalf));
                }

                if (domain.isInDomainCore(six) or domain.isInDomainCore(sevenandahalf)) {
                    CHECK(domain.isInDomainCore(six));
                    CHECK(domain.isInDomainCore(sevenandahalf));

                    CHECK(domain.isInDomainHalo(five));
                    CHECK_FALSE(domain.isInDomainHalo(twoandahalf));
                    CHECK_FALSE(domain.isInDomainCoreOrHalo(twoandahalf));
                    CHECK_FALSE(domain.isInDomainCore(twoandahalf));
                }
            }
            // check neighborhood

            if (rank != 0) {
                auto otherRank = rank == 1 ? 2 : 1;
                // neighborIndex (2,1,1) is the neighbor in the west (x+dx, y, z)
                CHECK(domain.neighborRanks()[domain.neighborIndex(1 + 1, 1, 1)] == otherRank);
                // (0, 1, 1) is the neighbor in the east
                CHECK(domain.neighborRanks()[domain.neighborIndex(1 - 1, 1, 1)] == otherRank);
                // "neighbor" (1,1,1) is this domain/rank
                CHECK(domain.neighborRanks()[domain.neighborIndex(1, 1, 1)] == rank);

                CHECK(domain.neighborTypes()[domain.neighborIndex(1 + 1, 1, 1)] == NeighborType::regular);
                CHECK(domain.neighborTypes()[domain.neighborIndex(1 - 1, 1, 1)] == NeighborType::regular);
                CHECK(domain.neighborTypes()[domain.neighborIndex(1, 1, 1)] == NeighborType::self);

                // all other neighbors do not exist -> ranks -1, types nan
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        for (int k = 0; k < 3; ++k) {
                            // skip the ones we checked above
                            if ((not(i == 1 and j == 1 and k == 1))
                                and (not(i == 0 and j == 1 and k == 1))
                                and (not(i == 2 and j == 1 and k == 1))) {
                                CHECK(domain.neighborRanks()[domain.neighborIndex(i, j, k)] == -1);
                                CHECK(domain.neighborTypes()[domain.neighborIndex(i, j, k)] == NeighborType::nan);
                            }
                        }
                    }
                }
            }
        }

        SECTION("Rong user given domain widths") {
            context.kernelConfiguration().mpi.dx = 2.2;
            context.kernelConfiguration().mpi.dy = 4.6;
            context.kernelConfiguration().mpi.dz = 4.6;
            context.kernelConfiguration().mpi.rank = 0;
            context.kernelConfiguration().mpi.worldSize = 4+1;
            auto ctor = [&]() {
                readdy::kernel::mpi::model::MPIDomain domain(context);
            };
            CHECK_THROWS(ctor());
        }
    }

    SECTION("Check one position in (4 by 2 by 1) domains, periodic in xy") {
        context.boxSize() = {20., 10., 1.};
        context.periodicBoundaryConditions() = {true, true, false};
        context.kernelConfiguration().mpi.dx = 4.6;
        context.kernelConfiguration().mpi.dy = 4.6;
        context.kernelConfiguration().mpi.dz = 4.6;
        int worldSize = 1 + (4 * 2 * 1);
        std::size_t posCount{0};
        readdy::Vec3 pos{-8., -4., 0.49};
        int south, north, west, east, northwest, southwest;
        for (int rank = 0; rank < worldSize; ++rank) {
            context.kernelConfiguration().mpi.rank = rank;
            context.kernelConfiguration().mpi.worldSize = worldSize;
            readdy::kernel::mpi::model::MPIDomain domain(context);
            CHECK(domain.worldSize() == 9);
            CHECK(domain.nDomainsPerAxis() == std::array<std::size_t, 3>({4, 2, 1}));
            CHECK(domain.domainIndex()(0, 0, 0) == 0);
            CHECK(domain.rank() == rank);
            CHECK(domain.haloThickness() == context.calculateMaxCutoff());

            for (const auto otherRank : domain.workerRanks()) {
                readdy::Vec3 origin, extent;
                std::tie(origin, extent) = domain.coreOfDomain(otherRank);
                auto center = origin + 0.5 * extent;
                CHECK(domain.rankOfPosition(center) == otherRank);
            }

            if (rank != 0) {
                if (domain.isInDomainCore(pos)) {
                    posCount++;

                    west = domain.neighborRanks()[domain.neighborIndex(0, 1, 1)];
                    east = domain.neighborRanks()[domain.neighborIndex(2, 1, 1)];
                    south = domain.neighborRanks()[domain.neighborIndex(1, 0, 1)];
                    north = domain.neighborRanks()[domain.neighborIndex(1, 2, 1)];
                    northwest = domain.neighborRanks()[domain.neighborIndex(0, 2, 1)];
                    southwest = domain.neighborRanks()[domain.neighborIndex(0, 0, 1)];

                    CHECK(northwest == southwest); // due to periodicity
                    CHECK(north == south); // due to periodicity
                }
            }
        }
        CHECK(posCount == 1); // can only be in one domain core

        context.kernelConfiguration().mpi.rank = south;
        readdy::kernel::mpi::model::MPIDomain southDomain(context);
        CHECK(southDomain.isInDomainHalo(pos));
        context.kernelConfiguration().mpi.rank = north;
        readdy::kernel::mpi::model::MPIDomain northDomain(context);
        CHECK(northDomain.isInDomainHalo(pos));
        context.kernelConfiguration().mpi.rank = west;
        readdy::kernel::mpi::model::MPIDomain westDomain(context);
        CHECK(westDomain.isInDomainHalo(pos));
        context.kernelConfiguration().mpi.rank = east;
        readdy::kernel::mpi::model::MPIDomain eastDomain(context);
        CHECK_FALSE(eastDomain.isInDomainHalo(pos));
        context.kernelConfiguration().mpi.rank = northwest;
        readdy::kernel::mpi::model::MPIDomain nwDomain(context);
        CHECK(nwDomain.isInDomainHalo(pos));
    }

    SECTION("Periodic in a direction (here y) which has only one domain") {
        context.boxSize() = {10., 5., 1.};
        context.periodicBoundaryConditions() = {true, true, false};
        // user given values should be at least twice the cutoff
        context.kernelConfiguration().mpi.dx = 4.6;
        context.kernelConfiguration().mpi.dy = 4.6;
        context.kernelConfiguration().mpi.dz = 4.6;
        int worldSize = 3;
        for (int rank = 0; rank < worldSize; ++rank) {
            context.kernelConfiguration().mpi.rank = rank;
            context.kernelConfiguration().mpi.worldSize = worldSize;
            readdy::kernel::mpi::model::MPIDomain domain(context);
            CHECK(domain.worldSize() == 3);
            CHECK(domain.nDomainsPerAxis() == std::array<std::size_t, 3>({2, 1, 1}));
            CHECK(domain.domainIndex()(0, 0, 0) == 0);
            CHECK(domain.rank() == rank);
            CHECK(domain.haloThickness() == context.calculateMaxCutoff());

            for (const auto otherRank : domain.workerRanks()) {
                readdy::Vec3 origin, extent;
                std::tie(origin, extent) = domain.coreOfDomain(otherRank);
                auto center = origin + 0.5 * extent;
                CHECK(domain.rankOfPosition(center) == otherRank);
            }

            if (rank != 0) {
                readdy::Vec3 five{4.99 - 5., -0.1, 0.};
                readdy::Vec3 twoandahalf{2.5 - 5., -0.1, 0.};
                readdy::Vec3 six{5.01 - 5., -0.1, 0.};
                readdy::Vec3 sevenandahalf{7.5 - 5., -0.1, 0.};

                if (domain.isInDomainCore(five) or domain.isInDomainCore(twoandahalf)) {
                    // if either of them is in domain, both must be in domain
                    CHECK(domain.isInDomainCore(five));
                    CHECK(domain.isInDomainCore(twoandahalf));

                    // six is not in this domain but very close
                    CHECK(domain.isInDomainHalo(six));
                    // sevenandahalf however is too far in the next domain, and thus not in this ones' halo, nor core
                    CHECK_FALSE(domain.isInDomainHalo(sevenandahalf));
                    CHECK_FALSE(domain.isInDomainCoreOrHalo(sevenandahalf));
                    CHECK_FALSE(domain.isInDomainCore(sevenandahalf));
                }

                if (domain.isInDomainCore(six) or domain.isInDomainCore(sevenandahalf)) {
                    CHECK(domain.isInDomainCore(six));
                    CHECK(domain.isInDomainCore(sevenandahalf));

                    CHECK(domain.isInDomainHalo(five));
                    CHECK_FALSE(domain.isInDomainHalo(twoandahalf));
                    CHECK_FALSE(domain.isInDomainCoreOrHalo(twoandahalf));
                    CHECK_FALSE(domain.isInDomainCore(twoandahalf));
                }
            }
            // check neighborhood
            if (rank != 0) {
                auto otherRank = rank == 1 ? 2 : 1;

                auto west = otherRank; // two boxes and periodic in WE direction
                auto east = otherRank;
                auto south = rank; // only one box and periodic in SN direction
                auto north = rank;

                CHECK(domain.neighborRanks()[domain.neighborIndex(1 - 1, 1, 1)] == west);
                CHECK(domain.neighborRanks()[domain.neighborIndex(1 + 1, 1, 1)] == east);
                CHECK(domain.neighborRanks()[domain.neighborIndex(1, 1 - 1, 1)] == south);
                CHECK(domain.neighborRanks()[domain.neighborIndex(1, 1 + 1, 1)] == north);

                // up and down are not periodic and very thin -> only one domain on this axis
                CHECK(domain.neighborRanks()[domain.neighborIndex(1, 1, 1 - 1)] == -1);
                CHECK(domain.neighborRanks()[domain.neighborIndex(1, 1, 1 + 1)] == -1);


                CHECK(domain.neighborRanks()[domain.neighborIndex(1, 1, 1)] == rank);

                CHECK(domain.neighborTypes()[domain.neighborIndex(1 + 1, 1, 1)] == NeighborType::regular);
                CHECK(domain.neighborTypes()[domain.neighborIndex(1 - 1, 1, 1)] == NeighborType::regular);
                CHECK(domain.neighborTypes()[domain.neighborIndex(1, 1 - 1, 1)] == NeighborType::self);
                CHECK(domain.neighborTypes()[domain.neighborIndex(1, 1 + 1, 1)] == NeighborType::self);
                CHECK(domain.neighborTypes()[domain.neighborIndex(1, 1, 1 - 1)] == NeighborType::nan);
                CHECK(domain.neighborTypes()[domain.neighborIndex(1, 1, 1 + 1)] == NeighborType::nan);
                CHECK(domain.neighborTypes()[domain.neighborIndex(1, 1, 1)] == NeighborType::self);

                // requesting a position outside the box makes no sense, except either along west or east
                auto center = domain.origin() + 0.5 * domain.extent();
                readdy::Vec3 dEW{3., 0., 0.};
                readdy::Vec3 dNS{0., 3., 0.};
                readdy::Vec3 dUD{0., 0., 3.};
                CHECK_THROWS(domain.rankOfPosition(center - dUD));
                CHECK_THROWS(domain.rankOfPosition(center + dUD));
                CHECK_THROWS(domain.rankOfPosition(center - dNS));
                CHECK_THROWS(domain.rankOfPosition(center + dNS));
                if (center[0] < 0.) {
                    CHECK_THROWS(domain.rankOfPosition(center - dEW));
                    CHECK(domain.rankOfPosition(center + dEW) == otherRank);
                } else {
                    CHECK_THROWS(domain.rankOfPosition(center + dEW));
                    CHECK(domain.rankOfPosition(center - dEW) == otherRank);
                }
            }
        }
    }

    SECTION("Test wrapIntoThisHalo") {
        std::vector<std::array<bool, 3>> pbcs{
                {false, false, false},
                {false, false, true},
                {false, true,  true},
                {true,  true,  true}};
        std::vector<std::pair<int, std::array<readdy::scalar, 3>>> boxes{
                {1 + 1,  {5.,  5.,  5.}},
                {2 + 1,  {10., 5.,  5.}},
                {4 + 1,  {10., 10., 5.}},
                {8 + 1,  {10., 10., 10.}},
                {16 + 1, {20., 10., 10.}},
                {32 + 1, {20., 20., 10.}},
                {12 + 1, {10., 15., 10.}}};

        context.kernelConfiguration().mpi.dx = 4.6;
        context.kernelConfiguration().mpi.dy = 4.6;
        context.kernelConfiguration().mpi.dz = 4.6;
        readdy::util::Index2D index2D(pbcs.size(), boxes.size());
        for (int i = 0; i < pbcs.size(); ++i) {
            const auto &pbc = pbcs[i];
            for (int j = 0; j < boxes.size(); ++j) {
                const auto &worldSize = boxes[j].first;
                const auto &box = boxes[j].second;
                DYNAMIC_SECTION("different pbcs and boxes " << index2D(i, j)) {
                    context.boxSize() = box;
                    context.periodicBoundaryConditions() = pbc;
                    // position is always at the positive "end" of the box
                    readdy::Vec3 pos{0.45 * box[0], 0.45 * box[1], 0.45 * box[2]};
                    for (int rank = 1; rank < worldSize; ++rank) {
                        context.kernelConfiguration().mpi.rank = rank;
                        context.kernelConfiguration().mpi.worldSize = worldSize;
                        readdy::kernel::mpi::model::MPIDomain domain(context);
                        auto wrapped = domain.wrapIntoThisHalo(pos);
                        if (domain.isInDomainCore(pos)) {
                            CHECK(pos == wrapped);
                        }
                        auto myIdx = domain.myIdx();
                        auto idx = domain.ijkOfPosition(pos);
                        for (int d = 0; d < 3; ++d) {
                            if (pbc[d]) {
                                if (idx[d] != myIdx[d]) {
                                    CHECK(pos[d] - wrapped[d] > 4.);
                                }
                            }
                        }

                        for (const auto otherRank : domain.workerRanks()) {
                            readdy::Vec3 origin, extent;
                            std::tie(origin, extent) = domain.coreOfDomain(otherRank);
                            auto center = origin + 0.5 * extent;
                            CHECK(domain.rankOfPosition(center) == otherRank);
                        }

                    }
                }
            }
        }
    }
}
