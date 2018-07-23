/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * << detailed description >>
 *
 * @file TestTimer.cpp
 * @brief << brief description >>
 * @author chrisfroe
 * @date 01.09.17
 * @copyright GPL-3
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