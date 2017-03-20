/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          * 
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


#include <array>

#include <gtest/gtest.h>

#include <readdy/common/tuple_utils.h>
#include <readdy/common/logging.h>
#include <readdy/common/ParticleTypeTuple.h>

/**
 * << detailed description >>
 *
 * @file TestTopologyGraphs.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

TEST(TestTopologyGraphs, TestQuadruple) {
    readdy::util::particle_type_quadruple_hasher hasher;
    std::array<readdy::particle_type_type, 4> range{1, 2, 3, 4};
    do {
        std::stringstream ss;
        ss << range[0] << ", " << range[1] << ", " << range[2] << ", " << range[3];
        EXPECT_EQ(readdy::util::sortTypeQuadruple(range[0], range[1], range[2], range[3]), std::make_tuple(1, 2, 3, 4))
                            << "failed for range " << ss.str() << ", should always yield a sorted tuple!";
    } while (std::next_permutation(range.begin(), range.end()));
    EXPECT_EQ(hasher(std::make_tuple(1, 2, 3, 4)), hasher(std::make_tuple(4, 3, 2, 1)));
    EXPECT_EQ(hasher(std::make_tuple(1, 3, 2, 4)), hasher(std::make_tuple(4, 2, 3, 1)));
}

TEST(TestTopologyGraphs, TestTriple) {
    readdy::util::particle_type_triple_hasher hasher;
    std::array<readdy::particle_type_type, 3> range{1, 2, 3};
    do {
        std::stringstream ss;
        ss << range[0] << ", " << range[1] << ", " << range[2];
        EXPECT_EQ(readdy::util::sortTypeTriple(range[0], range[1], range[2]), std::make_tuple(1, 2, 3))
                            << "failed for range " << ss.str() << ", should always yield a sorted tuple!";
    } while (std::next_permutation(range.begin(), range.end()));
    EXPECT_EQ(hasher(std::make_tuple(1, 2, 3)), hasher(std::make_tuple(3, 2, 1)));
    EXPECT_EQ(hasher(std::make_tuple(2, 1, 3)), hasher(std::make_tuple(3, 1, 2)));
}
