/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * @file TestAlgorithms.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 3/1/18
 */


#include <unordered_set>
#include "gtest/gtest.h"
#include "readdy/common/algorithm.h"


namespace {

using namespace readdy;

struct Event {
    int i;
    scalar rate;
    scalar cumulativeRate;
};

TEST(TestAlgorithms, EvaluateAll) {
    auto n = 1000U;
    std::vector<Event> events (n);
    for(auto i = 0U; i < n; ++i) {
        events.at(i).i = i;
        events.at(i).rate = i;
    }

    auto shouldEval = [](const Event &event) { return true; };
    auto depending = [](const Event &e1, const Event &e2) { return e1.rate == e2.rate; };

    std::unordered_set<int> set;
    auto eval = [&set](const Event &event) {
        set.insert(event.i);
    };

    algo::performEvents(events, shouldEval, depending, eval);

    for(auto i = 0; i < n; ++i) {
        ASSERT_TRUE(set.find(i) != set.end());
    }
}

TEST(TestAlgorithms, EvaluateEven) {
    auto n = 1000U;
    std::vector<Event> events (n);
    for(auto i = 0U; i < n; ++i) {
        events.at(i).i = i;
        events.at(i).rate = i/2;
    }

    auto shouldEval = [](const Event &event) { return true; };
    auto depending = [](const Event &e1, const Event &e2) { return e1.rate == e2.rate; };

    std::unordered_set<scalar> set;
    auto eval = [&set, &events](const Event &event) {
        {
            // check that the events are all moved to the end that were already taken care of
            for(const auto rate : set) {
                auto it1 = std::find_if(events.end() - 2*set.size(), events.end(), [rate](const Event &event) {
                    return event.rate == rate;
                });
                ASSERT_NE(it1, events.end());
                auto it2 = std::find_if(it1+1, events.end(), [rate](const Event &event) {
                    return event.rate == rate;
                });
                ASSERT_NE(it2, events.end());
            }
        }

        set.insert(event.rate);
    };

    algo::performEvents(events, shouldEval, depending, eval);

    for(auto i = 0; i < n/2; ++i) {
        ASSERT_TRUE(set.find(i) != set.end());
    }
}

}
