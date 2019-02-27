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
 * @copyright BSD-3
 */

#include <catch2/catch.hpp>

#include <readdy/common/index_persistent_vector.h>
#include <readdy/common/range.h>
#include <memory>

struct Element {
    int val;
    bool deactivated;

    Element() = default;
    Element(int val, bool deactivated) : val(val), deactivated(deactivated) {}

    void deactivate() {
        deactivated = true;
    }
};

struct NoElement {
    int val;
};

using ElementPtr = std::shared_ptr<Element>;

using namespace readdy;

TEMPLATE_TEST_CASE("Test the index persistent vector class.", "[index-persistent-vector]", Element, ElementPtr) {
    util::index_persistent_vector<Element> vec;

    SECTION("Can be deactivated") {
        CHECK(util::detail::can_be_deactivated<Element>::value);
        CHECK_FALSE(util::detail::can_be_deactivated<NoElement>::value);
        CHECK_FALSE(util::detail::can_be_deactivated<std::string>::value);
    }

    SECTION("Add and remove") {
        REQUIRE(vec.empty());
        vec.push_back({});
        vec.push_back({});
        vec.push_back({});

        REQUIRE(vec.size() == 3);
        REQUIRE_FALSE(vec.empty());
        vec.erase(vec.begin());
        vec.erase(vec.begin()+1);
        vec.erase(vec.begin()+2);

        REQUIRE(vec.size() == 3);
        REQUIRE(vec.n_deactivated() == 3);
        REQUIRE(vec.empty());

        for(const auto& x : vec) {
            CHECK(x.deactivated);
        }
    }

    SECTION("Erase range") {
        util::range<int> range {0, 10};
        std::for_each(range.begin(), range.end(), [&](int) { vec.push_back({}); });

        REQUIRE(vec.size() == 10);
        vec.erase(vec.begin(), vec.end()-5);

        REQUIRE(vec.size() == 10);
        REQUIRE(vec.n_deactivated() == 5);

        for(std::size_t i = 0; i < vec.size(); ++i) {
            if(i < 5) {
                REQUIRE(vec.at(i).deactivated);
            } else {
                REQUIRE_FALSE(vec.at(i).deactivated);
            }
        }
    }

    SECTION("Reclaim index") {
        vec.push_back({});
        vec.push_back({});

        vec.erase(vec.begin());
        REQUIRE(vec.size() == 2);
        REQUIRE(vec.n_deactivated() == 1);

        vec.emplace_back(-1, false);
        REQUIRE(vec.size() == 2);
        REQUIRE(vec.n_deactivated() == 0);
        REQUIRE(vec.begin()->val == -1);
        REQUIRE_FALSE(vec.begin()->deactivated);

        for(const auto& x : vec) {
            CHECK_FALSE(x.deactivated);
        }
    }

    SECTION("Deactivation of elements") {
        vec.push_back({});
        REQUIRE(vec.size() == 1);
        REQUIRE(vec.n_deactivated() == 0);
        REQUIRE_FALSE(vec.begin()->deactivated);

        vec.erase(vec.begin());
        REQUIRE(vec.size() == 1);
        REQUIRE(vec.n_deactivated() == 1);
        REQUIRE(vec.begin()->deactivated);

        vec.push_back({});
        REQUIRE(vec.size() == 1);
        REQUIRE(vec.n_deactivated() == 0);
        REQUIRE_FALSE(vec.begin()->deactivated);
    }
}
