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

#include <gtest/gtest.h>
#include <readdy/common/common.h>

namespace {

using vec = readdy::Vec3;

TEST(Vec3, SizeOfVec3) {
    vec v(0, 0, 0);
    EXPECT_EQ(sizeof(readdy::scalar) * 3, sizeof(v)) << "a vector should have exactly 24 bytes";
}

TEST(Vec3, PlusEq) {
    vec v(1, 2, 3);
    vec v2(1, 1, 1);
    auto v3 = v + v2;
    v += v2;
    EXPECT_EQ(v, vec(2, 3, 4));
    EXPECT_EQ(v, v3);
}

TEST(Vec3, MinusEq) {
    vec v(1.1, 2, 3);
    vec v2(1, 1, 1);
    auto v3 = v - v2;
    v -= v2;
    EXPECT_TRUE(v.almostEquals(vec(0.1, 1, 2)));
    EXPECT_EQ(v, v3);
}

TEST(Vec3, PlusEqScalar) {
    vec v(1, 2, 3);
    v += 1.;
    EXPECT_EQ(v, vec(2, 3, 4));
}

TEST(Vec3, MinusEqScalar) {
    vec v(1, 2, 3);
    v -= 1.;
    EXPECT_EQ(v, vec(0, 1, 2));
}

TEST(Vec3, Scale) {
    vec v(1, 2, 3);
    v *= .5;
    EXPECT_EQ(v, vec(.5, 1, 1.5));
}

TEST(Vec3, Divide) {
    vec v(1, 2, 3);
    auto v2 = v / 2;
    v /= 2.;
    EXPECT_EQ(v, vec(.5, 1., 1.5));
    EXPECT_EQ(v, v2);
}

TEST(Vec3, ModifyInplace) {
    vec v(1, 2, 3);
    vec v2(1, 2, 3);
    vec v3(1, 2, 3);
    v += -1 * v3;
    v2 -= v3;
    EXPECT_EQ(v, vec(0, 0, 0));
    EXPECT_EQ(v, v2);
    EXPECT_EQ(v3, vec(1, 2, 3));
}

TEST(Vec3, GEQ) {
    vec v(1, 2, 3);
    vec v2(1, 2, 3);
    vec v3(2, 3, 4);
    EXPECT_TRUE(v >= v2);
    EXPECT_FALSE(v > v2);
    EXPECT_TRUE(v3 > v);
    EXPECT_FALSE(v3 < v);
    EXPECT_GE(v, v2);
    EXPECT_GT(v3, v);
}

TEST(Vec3, LEQ) {
    vec v(1, 2, 3);
    vec v2(1, 2, 3);
    vec v3(0, 1, 2);
    EXPECT_TRUE(v <= v2);
    EXPECT_FALSE(v < v2);
    EXPECT_TRUE(v3 < v);
    EXPECT_FALSE(v < v3);
    EXPECT_LE(v, v2);
    EXPECT_LT(v3, v);
}

TEST(Vec3, InnerProduct) {
    vec v(1, 2, 3);
    vec v2(2, 3, 4);
    EXPECT_EQ(v * v2, 1 * 2 + 2 * 3 + 3 * 4);
}

TEST(Vec3, ProductWithScalar) {
    vec v(1, 2, 3);
    vec v2 = v;
    v2 *= 5;
    EXPECT_EQ(5 * v, vec(5, 10, 15));
    EXPECT_EQ(v * 5, vec(5, 10, 15));
    EXPECT_EQ(v2, 5 * v);
}

TEST(Vec3, DivideByScalar) {
    vec v(3, 4, 5);
    vec v2 = v;
    v2 /= 2;
    EXPECT_EQ(v / 2, vec(1.5, 2, 2.5));
    EXPECT_EQ(v2, v / 2);
}

TEST(Vec3, MinusScalar) {
    vec v(1, 2, 3);
    vec v2 = v;
    v2 -= 1;
    EXPECT_EQ(v - 1, vec(0, 1, 2));
    EXPECT_EQ(v - 1, v2);
}

TEST(Vec3, Norm) {
    vec v(2, 2, 2);
    EXPECT_EQ(v.norm(), std::sqrt(2 * 2 + 2 * 2 + 2 * 2));
}

TEST(Vec3, NormSquared) {
    vec v(-.9, .1, .1);
    EXPECT_EQ(v.normSquared(), .9 * .9 + .1 * .1 + .1 * .1);
    EXPECT_EQ(v.normSquared(), v*v);
}

}
