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
 * @file TestVec3.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#include <catch2/catch_test_macros.hpp>
#include <readdy/common/common.h>

using vec = readdy::Vec3;

TEST_CASE("Test vec3", "[vec3]") {
    SECTION("Size of vec3") {
        vec v(0, 0, 0);
        // vec3 should be behaving like POD
        REQUIRE(sizeof(readdy::scalar) * 3 == sizeof(v));
    }
    SECTION("+=") {
        vec v(1, 2, 3);
        vec v2(1, 1, 1);
        auto v3 = v + v2;
        v += v2;
        REQUIRE(v == vec(2, 3, 4));
        REQUIRE(v == v3);
    }
    SECTION("-=") {
        vec v(1.1, 2, 3);
        vec v2(1, 1, 1);
        auto v3 = v - v2;
        v -= v2;
        REQUIRE(v.almostEquals(vec(0.1, 1, 2)));
        REQUIRE(v == v3);
    }
    SECTION("+= scalar") {
        vec v(1, 2, 3);
        v += 1.;
        REQUIRE(v == vec(2, 3, 4));
    }
    SECTION("-= scalar") {
        vec v(1, 2, 3);
        v -= 1.;
        REQUIRE(v == vec(0, 1, 2));
    }
    SECTION("scale") {
        vec v(1, 2, 3);
        v *= .5;
        REQUIRE(v == vec(.5, 1, 1.5));
    }
    SECTION("divide") {
        vec v(1, 2, 3);
        auto v2 = v / 2;
        v /= 2.;
        REQUIRE(v == vec(.5, 1., 1.5));
        REQUIRE(v == v2);
    }
    SECTION("In-place modification") {
        vec v(1, 2, 3);
        vec v2(1, 2, 3);
        vec v3(1, 2, 3);
        v += -1 * v3;
        v2 -= v3;
        REQUIRE(v == vec(0, 0, 0));
        REQUIRE(v == v2);
        REQUIRE(v3 == vec(1, 2, 3));
    }
    SECTION(">=") {
        vec v(1, 2, 3);
        vec v2(1, 2, 3);
        vec v3(2, 3, 4);
        REQUIRE(v >= v2);
        REQUIRE_FALSE(v > v2);
        REQUIRE(v3 > v);
        REQUIRE_FALSE(v3 < v);
        REQUIRE(v >= v2);
        REQUIRE(v3 > v);
    }
    SECTION("<=") {
        vec v(1, 2, 3);
        vec v2(1, 2, 3);
        vec v3(0, 1, 2);
        REQUIRE(v <= v2);
        REQUIRE_FALSE(v < v2);
        REQUIRE(v3 < v);
        REQUIRE_FALSE(v < v3);
        REQUIRE(v <= v2);
        REQUIRE(v3 < v);
    }
    SECTION("< . , . >") {
        vec v(1, 2, 3);
        vec v2(2, 3, 4);
        REQUIRE(v * v2 == 1 * 2 + 2 * 3 + 3 * 4);
    }
    SECTION("-= scalar") {
        vec v(1, 2, 3);
        vec v2 = v;
        v2 -= 1;
        REQUIRE(v - 1 == vec(0, 1, 2));
        REQUIRE(v - 1 == v2);
    }
    SECTION("|| . ||") {
        vec v(2, 2, 2);
        REQUIRE(v.norm() == std::sqrt(2 * 2 + 2 * 2 + 2 * 2));
    }
    SECTION("|| . ||^2") {
        vec v(-.9, .1, .1);
        REQUIRE(v.normSquared() == .9 * .9 + .1 * .1 + .1 * .1);
        REQUIRE(v.normSquared() == v*v);
    }
}
