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
