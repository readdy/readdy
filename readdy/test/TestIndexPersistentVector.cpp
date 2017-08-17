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


/**
 * << detailed description >>
 *
 * @file TestIndexPersistentVector.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/common/index_persistent_vector.h>
#include <readdy/common/range.h>
#include <memory>
#include "gtest/gtest.h"

namespace {

struct Element {
    int val;
    bool deactivated;

    Element() = default;
    Element(int val, bool deactivated) : val(val), deactivated(deactivated) {}

    void deactivate() {
        deactivated = true;
    }
};

TEST(TestIndexPersistentVector, DetectDeactivateable) {
    using namespace readdy;

    struct NoElement {
        int val;
    };

    ASSERT_TRUE(util::detail::can_be_deactivated<Element>::value);
    ASSERT_FALSE(util::detail::can_be_deactivated<NoElement>::value);
    ASSERT_FALSE(util::detail::can_be_deactivated<std::string>::value);
}

TEST(TestIndexPersistentVector, AddAndRemove) {
    using namespace readdy;
    util::index_persistent_vector<Element> vec;

    ASSERT_TRUE(vec.empty());
    vec.push_back({});
    vec.push_back({});
    vec.push_back({});

    ASSERT_EQ(vec.size(), 3);
    ASSERT_FALSE(vec.empty());
    vec.erase(vec.begin());
    vec.erase(vec.begin()+1);
    vec.erase(vec.begin()+2);

    ASSERT_EQ(vec.size(), 3);
    ASSERT_EQ(vec.n_deactivated(), 3);
    ASSERT_TRUE(vec.empty());

    for(const auto& x : vec) {
        ASSERT_TRUE(x.deactivated);
    }
}

TEST(TestIndexPersistentVector, EraseRange) {
    using namespace readdy;
    util::index_persistent_vector<Element> vec;
    util::range<int> range {0, 10};
    std::for_each(range.begin(), range.end(), [&](int) { vec.push_back({}); });

    ASSERT_EQ(vec.size(), 10);
    vec.erase(vec.begin(), vec.end()-5);

    ASSERT_EQ(vec.size(), 10);
    ASSERT_EQ(vec.n_deactivated(), 5);

    for(std::size_t i = 0; i < vec.size(); ++i) {
        if(i < 5) {
            EXPECT_TRUE(vec.at(i).deactivated);
        } else {
            EXPECT_FALSE(vec.at(i).deactivated);
        }
    }
}

TEST(TestIndexPersistentVector, ReclaimIndex) {
    using namespace readdy;
    util::index_persistent_vector<Element> vec;
    vec.push_back({});
    vec.push_back({});

    vec.erase(vec.begin());
    ASSERT_EQ(vec.size(), 2);
    ASSERT_EQ(vec.n_deactivated(), 1);

    vec.emplace_back(-1, false);
    ASSERT_EQ(vec.size(), 2);
    ASSERT_EQ(vec.n_deactivated(), 0);
    ASSERT_EQ(vec.begin()->val, -1);
    ASSERT_EQ(vec.begin()->deactivated, false);

    for(const auto& x : vec) {
        ASSERT_FALSE(x.deactivated);
    }
}

TEST(TestIndexPersistentVector, Deactivation) {
    using namespace readdy;
    util::index_persistent_vector<Element> vec;

    vec.push_back({});
    ASSERT_EQ(vec.size(), 1);
    ASSERT_EQ(vec.n_deactivated(), 0);
    ASSERT_FALSE(vec.begin()->deactivated);

    vec.erase(vec.begin());
    ASSERT_EQ(vec.size(), 1);
    ASSERT_EQ(vec.n_deactivated(), 1);
    ASSERT_TRUE(vec.begin()->deactivated);

    vec.push_back({});
    ASSERT_EQ(vec.size(), 1);
    ASSERT_EQ(vec.n_deactivated(), 0);
    ASSERT_FALSE(vec.begin()->deactivated);
}

TEST(TestIndexPersistentVector, PointerDeactivation) {
    using namespace readdy;
    util::index_persistent_vector<std::shared_ptr<Element>> vec;

    vec.push_back(std::make_shared<Element>());
    ASSERT_EQ(vec.size(), 1);
    ASSERT_EQ(vec.n_deactivated(), 0);
    ASSERT_FALSE((*vec.begin())->deactivated);

    vec.erase(vec.begin());
    ASSERT_EQ(vec.size(), 1);
    ASSERT_EQ(vec.n_deactivated(), 1);
    ASSERT_TRUE((*vec.begin())->deactivated);

    vec.push_back(std::make_shared<Element>());
    ASSERT_EQ(vec.size(), 1);
    ASSERT_EQ(vec.n_deactivated(), 0);
    ASSERT_FALSE((*vec.begin())->deactivated);
}

}