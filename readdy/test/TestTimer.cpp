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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the     *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file TestTimer.cpp
 * @brief << brief description >>
 * @author chrisfroe
 * @date 01.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <gtest/gtest.h>
#include <readdy/common/Timer.h>

namespace {

struct TestTimer : ::testing::Test {
};

using data = readdy::util::PerformanceData;
using node = readdy::util::PerformanceNode;

TEST(TestTimer, TimerMeasureOnce) {
    node n("label", true);
    {
        auto t = n.timeit();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    EXPECT_GT(n.data().cumulativeTime, static_cast<data::time>(0));
    EXPECT_EQ(n.data().count, 1);
}

TEST(TestTimer, Clear) {
    node n1("1", true);
    node n2("2", true);
    {
        auto t1 = n1.timeit();
        auto t2 = n2.timeit();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    n2.clear();
    {
        auto t1 = n1.timeit();
        auto t2 = n2.timeit();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    EXPECT_EQ(n1.data().count, 2);
    EXPECT_EQ(n2.data().count, 1) << "n2 was cleared in between, thus only 1 call";
}

}