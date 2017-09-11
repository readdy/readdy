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
 * @file TextIndex.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include "gtest/gtest.h"
#include <readdy/common/Index.h>

namespace {

TEST(TestIndex, Test1D) {
    using namespace readdy::util;

    Index1D index (5);
    ASSERT_EQ(index[0], 5);
    ASSERT_EQ(index.get<0>(), 5);
    ASSERT_EQ(index.nElements(), 5);
    ASSERT_EQ(index(3), 3);
    ASSERT_EQ(index.inverse(3)[0], 3);
}
TEST(TestIndex, Test2D) {
    using namespace readdy::util;
    auto nrows = 5;
    auto ncols = 6;
    Index2D index (nrows, ncols);
    ASSERT_EQ(index[0], nrows);
    ASSERT_EQ(index[1], ncols);
    ASSERT_EQ(index.get<0>(), nrows);
    ASSERT_EQ(index.get<1>(), ncols);
    ASSERT_EQ(index.nElements(), nrows * ncols);
    ASSERT_EQ(index(0, 3), 3);
    ASSERT_EQ(index(1, 3), 3 + ncols);
    ASSERT_EQ(index(4, 5), 5 + 4 * ncols);
    ASSERT_EQ(index.inverse(3)[0], 0);
    ASSERT_EQ(index.inverse(3)[1], 3);
    ASSERT_EQ(index.inverse(3+ncols)[0], 1);
    ASSERT_EQ(index.inverse(3+ncols)[1], 3);
    ASSERT_EQ(index.inverse(5+4*ncols)[0], 4);
    ASSERT_EQ(index.inverse(5+4*ncols)[1], 5);
}
TEST(TestIndex, Test3D) {
    using namespace readdy::util;
    auto width = 6;
    auto height = 5;
    auto depth = 6;
    Index3D index (width, height, depth);
    ASSERT_EQ(index[0], index.get<0>());
    ASSERT_EQ(index[1], index.get<1>());
    ASSERT_EQ(index[2], index.get<2>());
    ASSERT_EQ(index[0], width);
    ASSERT_EQ(index[1], height);
    ASSERT_EQ(index[2], depth);
    ASSERT_EQ(index.nElements(), width * height * depth);
    ASSERT_EQ(index(0, 0, 3), 3);
    ASSERT_EQ(index(1, 2, 3), 3 + depth * (2 + height * 1));
    Index3D::GridDims expected {{1,2,3}};
    ASSERT_EQ(index.inverse(3 + depth * (2 + height * 1)), expected);
}

}
