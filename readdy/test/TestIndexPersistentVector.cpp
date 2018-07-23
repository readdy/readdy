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
 * @file TestIndexPersistentVector.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.17
 * @copyright GPL-3
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