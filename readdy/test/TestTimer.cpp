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

using timer = readdy::util::Timer;
using raii_timer = readdy::util::RAIITimer;

TEST(TestTimer, RAIITimerMeasureOnce) {
    {
        raii_timer t(true, "label");
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    EXPECT_GT(timer::times().at("label"), static_cast<timer::time>(0));
    EXPECT_EQ(timer::counts().at("label"), 1);
}

TEST(TestTimer, RAIITimerMeasureMultiple) {
    {
        raii_timer t1(true, "label1");
        raii_timer t2(true, "label2");
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    {
        raii_timer t3(true, "label2");
        raii_timer t4(true, "label2");
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    EXPECT_GT(timer::times().at("label1"), static_cast<timer::time>(0));
    EXPECT_EQ(timer::counts().at("label1"), 1);
    EXPECT_GT(timer::times().at("label2"), static_cast<timer::time>(0));
    EXPECT_EQ(timer::counts().at("label2"), 3);
}

TEST(TestTimer, Clear) {
    {
        raii_timer t(true, "label");
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    timer::clear();
    {
        raii_timer t(true, "label");
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    EXPECT_EQ(timer::counts().at("label"), 1);
}

}