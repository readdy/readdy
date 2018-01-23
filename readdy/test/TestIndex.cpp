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
#include <readdy/common/common.h>
#include <readdy/common/Index.h>

namespace {

TEST(TestIndex, Test1D) {
    using namespace readdy;

    util::Index1D index (5_z);
    ASSERT_EQ(index[0], 5);
    ASSERT_EQ(index.get<0>(), 5);
    ASSERT_EQ(index.nElements(), 5);
    ASSERT_EQ(index(3_z), 3);
    ASSERT_EQ(index.inverse(3)[0], 3);
}
TEST(TestIndex, Test2D) {
    using namespace readdy;
    auto nrows = 5_z;
    auto ncols = 6_z;
    util::Index2D index (nrows, ncols);
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
    using namespace readdy;
    auto width = 6_z;
    auto height = 5_z;
    auto depth = 6_z;
    util::Index3D index (width, height, depth);
    ASSERT_EQ(index[0], index.get<0>());
    ASSERT_EQ(index[1], index.get<1>());
    ASSERT_EQ(index[2], index.get<2>());
    ASSERT_EQ(index[0], width);
    ASSERT_EQ(index[1], height);
    ASSERT_EQ(index[2], depth);
    ASSERT_EQ(index.nElements(), width * height * depth);
    ASSERT_EQ(index(0, 0, 3), 3);
    ASSERT_EQ(index(1, 2, 3), 3 + depth * (2 + height * 1));
    util::Index3D::GridDims expected {{1_z,2_z,3_z}};
    ASSERT_EQ(index.inverse(3 + depth * (2 + height * 1)), expected);
}

TEST(TestIndex, Test3DElch) {
    using namespace readdy;
    auto width = 7_z;
    auto height = 13_z;
    auto depth = 5_z;

    util::Index3D index (width, height, depth);

    int n = 0;
    for(int i = 0; i < width; ++i) {
        for(int j = 0; j < height; ++j) {
            for(int k = 0; k < depth; ++k) {
                auto x = index(i, j, k);
                ASSERT_EQ(x, n);
                util::Index3D::GridDims expectedInverse {
                        {static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(k)}
                };
                ASSERT_EQ(index.inverse(x), expectedInverse);
                ++n;
            }
        }
    }
}

TEST(TestIndex, Test5DElch) {
    using namespace readdy;
    auto d1 = 7_z;
    auto d2 = 13_z;
    auto d3 = 5_z;
    auto d4 = 11_z;
    auto d5 = 23_z;

    util::Index<5> index (d1, d2, d3, d4, d5);

    int n = 0;
    for(int i1 = 0; i1 < d1; ++i1) {
        for(int i2 = 0; i2 < d2; ++i2) {
            for(int i3 = 0; i3 < d3; ++i3) {
                for(int i4 = 0; i4 < d4; ++i4) {
                    for(int i5 = 0; i5 < d5; ++i5) {
                        auto x = index(i1, i2, i3, i4, i5);
                        ASSERT_EQ(x, n);
                        util::Index<5>::GridDims expectedInverse {
                                {
                                        static_cast<std::size_t>(i1), static_cast<std::size_t>(i2),
                                        static_cast<std::size_t>(i3), static_cast<std::size_t>(i4),
                                        static_cast<std::size_t>(i5)
                                }
                        };
                        ASSERT_EQ(index.inverse(x), expectedInverse);
                        ++n;
                    }
                }
            }
        }
    }
}

}
