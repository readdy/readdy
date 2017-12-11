/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file TestNeighborListIterator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 15.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include "gtest/gtest.h"

#include "readdy/common/common.h"
#include "readdy/kernel/cpu_legacy/nl/NeighborListIterator.h"

namespace {

TEST(TestNeighborListIterator, Adaptive) {
    using namespace readdy;
    std::vector<std::vector<std::size_t>> adaptive_data;
    adaptive_data.emplace_back();
    adaptive_data.emplace_back(std::vector<std::size_t>{1_z, 2_z, 3_z});
    adaptive_data.emplace_back(std::vector<std::size_t>{3_z, 4_z});
    adaptive_data.emplace_back();
    adaptive_data.emplace_back(std::vector<std::size_t>{5_z});

    kernel::cpu::nl::NeighborListIterator itBegins(adaptive_data.begin(), adaptive_data.end(), true);
    kernel::cpu::nl::NeighborListIterator itEnds (adaptive_data.end(), adaptive_data.end(), true);

    auto it_proper = adaptive_data.begin();

    for(auto it = itBegins; it != itEnds; ++it, ++it_proper) {
        auto pit_proper = it_proper->begin();
        for(auto pit = it->begin(); pit != it->end(); ++pit, ++pit_proper) {
            EXPECT_EQ(*pit_proper, *pit);
        }
    }
}

TEST(TestNeighborListIterator, EmptyStatic) {
    using namespace readdy;
    std::vector<std::vector<std::size_t>> static_data;

    auto itBegins = kernel::cpu::nl::NeighborListIterator(static_data.begin(), static_data.end(), false);
    auto itEnds = kernel::cpu::nl::NeighborListIterator(static_data.end(), static_data.end(), false);

    ASSERT_EQ(itBegins, itEnds);
}

TEST(TestNeighborListIterator, EmptyDynamic) {
    using namespace readdy;
    std::vector<std::vector<std::size_t>> dynamic_data;

    auto itBegins = kernel::cpu::nl::NeighborListIterator(dynamic_data.begin(), dynamic_data.end(), true);
    auto itEnds = kernel::cpu::nl::NeighborListIterator(dynamic_data.end(), dynamic_data.end(), true);

    ASSERT_EQ(itBegins, itEnds);
}

TEST(TestNeighborListIterator, AlmostEmptyStatic) {
    using namespace readdy;
    std::vector<std::vector<std::size_t>> static_data {{}, {}, {}, {0_z, 0_z, 1_z, 3_z, 4_z, 0_z, 2_z}, {}, {}};

    auto itBegins = kernel::cpu::nl::NeighborListIterator(static_data.begin(), static_data.end(), false);
    auto itEnds = kernel::cpu::nl::NeighborListIterator(static_data.end(), static_data.end(), false);

    for(auto it = itBegins; it != itEnds; ++it) {
        EXPECT_TRUE(it->current_particle() == 0 || it->current_particle() == 1);
        if(it->current_particle() == 0) {
            EXPECT_EQ(it->n_neighbors(), 0);
        } else {
            EXPECT_EQ(it->n_neighbors(), 3);
            auto n1 = *it->begin();
            auto n2 = *(it->begin() + 1);
            auto n3 = *(it->begin() + 2);
            EXPECT_EQ(n1, 4);
            EXPECT_EQ(n2, 0);
            EXPECT_EQ(n3, 2);
        }
    }
}

TEST(TestNeighborListIterator, Static) {
    using namespace readdy;
    std::vector<std::vector<std::size_t>> static_data;
    // particle 0: no neighbors, particle 1: neighbors (2), particle 2: neighbors(4, 0, 1)
    static_data.emplace_back(std::vector<std::size_t>{0_z, 0_z, 1_z, 1_z, 2_z, 2_z, 3_z, 4_z, 0_z, 1_z});
    // this thread was lazy
    static_data.emplace_back(std::vector<std::size_t>{});
    // particle 3: neighbors (4), particle 4: no neighbors, particle 5: neighbors (1, 0)
    static_data.emplace_back(std::vector<std::size_t>{4_z, 0_z, 5_z, 2_z, 1_z, 0_z, 3_z, 1_z, 4_z});
    // this thread was lazy again
    static_data.emplace_back();
    static_data.emplace_back();

    kernel::cpu::nl::NeighborListIterator itBegins (static_data.begin(), static_data.end(), false);
    kernel::cpu::nl::NeighborListIterator itEnds (static_data.end(), static_data.end(), false);

    std::vector<std::size_t> expectedNNeighbors {0_z, 1_z, 3_z, 1_z, 0_z, 2_z};
    std::vector<std::vector<std::size_t>> expectedNeighbors {
            std::vector<std::size_t> {},
            std::vector<std::size_t> {2_z},
            std::vector<std::size_t> {4_z, 0_z, 1_z},
            std::vector<std::size_t> {4_z},
            std::vector<std::size_t> {},
            std::vector<std::size_t> {1_z, 0_z}
    };

    for(auto it = itBegins; it != itEnds; ++it) {
        //log::warn("Particle {} got {} neighbors", it->current_particle(), it->n_neighbors());
        EXPECT_EQ(it->n_neighbors(), expectedNNeighbors.at(it->current_particle()));
        EXPECT_EQ(expectedNeighbors.at(it->current_particle()).size(), it->n_neighbors());
        auto itExpected = expectedNeighbors.at(it->current_particle()).begin();
        for(auto neighbor : *it) {
            EXPECT_EQ(neighbor, *itExpected);
            ++itExpected;
        }
    }

    {
        auto itBeginsDynamic = kernel::cpu::nl::NeighborListIterator(expectedNeighbors.begin(), expectedNeighbors.end(), true);
        auto itEndsDynamic = kernel::cpu::nl::NeighborListIterator(expectedNeighbors.end(), expectedNeighbors.end(), true);

        auto itStatic = itBegins;
        const auto &itDynamic = itBeginsDynamic;
        for(; itStatic != itEnds; ++itStatic) {
            auto itd = itDynamic + itStatic->current_particle();
            EXPECT_EQ(itStatic->n_neighbors(), itd->n_neighbors());
            auto it1 = itStatic->begin();
            auto it2 = itd->begin();
            for(; it1 != itStatic->end(); ++it1, ++it2) {
                EXPECT_EQ(*it1, *it2);
            }
        }
    }
}

}