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
#include "readdy/kernel/cpu/nl/NeighborListIterator.h"

namespace {

TEST(TestNeighborListIterator, Adaptive) {
    using namespace readdy;
    std::vector<std::vector<std::size_t>> adaptive_data;
    adaptive_data.emplace_back();
    adaptive_data.emplace_back(std::vector<std::size_t>{1_z, 2_z, 3_z});
    adaptive_data.emplace_back(std::vector<std::size_t>{3_z, 4_z});
    adaptive_data.emplace_back();
    adaptive_data.emplace_back(std::vector<std::size_t>{5_z});

    kernel::cpu::nl::NeighborListIterator itBegins(adaptive_data.begin());
    kernel::cpu::nl::NeighborListIterator itEnds (adaptive_data.end());

    auto it_proper = adaptive_data.begin();

    for(auto it = itBegins; it != itEnds; ++it, ++it_proper) {
        auto pit_proper = it_proper->begin();
        for(auto pit = it->begin(); pit != it->end(); ++pit, ++pit_proper) {
            EXPECT_EQ(*pit_proper, *pit);
        }
    }
}

TEST(TestNeighborListIterator, Static) {

}

}