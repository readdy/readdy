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
#include <readdy/common/thread/executor.h>
#include <readdy/common/thread/Config.h>

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
    EXPECT_GT(n.data().cumulativeTime(), static_cast<data::time>(0));
    EXPECT_EQ(n.data().count(), 1);
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
    EXPECT_EQ(n1.data().count(), 2);
    EXPECT_EQ(n2.data().count(), 1) << "n2 was cleared in between, thus only one call";
}

TEST(TestTimer, Threaded) {
    node n("knoten", true);
    auto worker = [](std::size_t, node& nn){
        auto t = nn.timeit();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    };
    {
        readdy::util::thread::Config config;
        config.setNThreads(3);
        config.setMode(readdy::util::thread::ThreadMode::std_thread);
        const auto& executor = *config.executor();
        std::vector<std::function<void(std::size_t)>> executables;
        executables.reserve(config.nThreads());
        for (auto i = 0; i < config.nThreads(); ++i) {
            executables.push_back(executor.pack(worker, std::ref(n)));
        }
        executor.execute_and_wait(std::move(executables));
    }
    EXPECT_EQ(n.data().count(), 3);
    EXPECT_GT(n.data().cumulativeTime(), static_cast<data::time>(0));
}

TEST(TestTimer, SlashedPath) {
    node top("top", false);
    auto &mid = top.subnode("mid");
    auto &bot = mid.subnode("bot");
    EXPECT_EQ(top.child("mid /bot ").name(), "bot");
}

TEST(TestTimer, InvalidNodeName) {
    auto createNode = [](){
        node n("lk/afasov", false);
    };
    auto createSubNode = [](){
        node n("validname", false);
        auto &nn = n.subnode("invalid/name");
    };
    auto createNodeWhitespace = [](){
        node n(" hasleadingwhitespace", false);
    };
    EXPECT_THROW(createNode(), std::invalid_argument);
    EXPECT_THROW(createSubNode(), std::invalid_argument);
    EXPECT_THROW(createNodeWhitespace(), std::invalid_argument);
}

TEST(TestTimer, GetRootFromChild) {
    node top("top", false);
    auto &bottom = top.subnode("bottom");
    EXPECT_EQ(bottom.child("/").name(), "top");
}

}