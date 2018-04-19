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
 * @file TestMatrix33.cpp
 * @author clonker
 * @date 1/17/18
 */

#include "gtest/gtest.h"
#include <readdy/common/numeric.h>
#include <readdy/common/common.h>

using mat = readdy::Matrix33;
using vec = readdy::Vec3;

namespace {

TEST(Matrix33, OuterProduct) {
    vec u(1, 2, 3);
    vec v(3, 2, 1);
    auto outer = readdy::math::outerProduct<mat>(u, v);
    EXPECT_EQ(outer.at(0, 0), u[0] * v[0]);
    EXPECT_EQ(outer.at(0, 1), u[0] * v[1]);
    EXPECT_EQ(outer.at(0, 2), u[0] * v[2]);
    EXPECT_EQ(outer.at(1, 0), u[1] * v[0]);
    EXPECT_EQ(outer.at(1, 1), u[1] * v[1]);
    EXPECT_EQ(outer.at(1, 2), u[1] * v[2]);
    EXPECT_EQ(outer.at(2, 0), u[2] * v[0]);
    EXPECT_EQ(outer.at(2, 1), u[2] * v[1]);
    EXPECT_EQ(outer.at(2, 2), u[2] * v[2]);
}

TEST(Matrix33, Access) {
    mat m{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
    EXPECT_EQ(m.at(0, 0), 1);
    EXPECT_EQ(m.at(0, 1), 2);
    EXPECT_EQ(m.at(0, 2), 3);
    EXPECT_EQ(m.at(1, 0), 4);
    EXPECT_EQ(m.at(1, 1), 5);
    EXPECT_EQ(m.at(1, 2), 6);
    EXPECT_EQ(m.at(2, 0), 7);
    EXPECT_EQ(m.at(2, 1), 8);
    EXPECT_EQ(m.at(2, 2), 9);
}

TEST(Matrix33, Dims) {
    EXPECT_EQ(mat::n(), 3);
    EXPECT_EQ(mat::m(), 3);
}

TEST(Matrix33, Copy) {
    mat m{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
    auto copy = m;
    mat copy2(m);

    copy.at(0, 1) = 5;
    copy2.at(0, 1) = 6;

    EXPECT_EQ(m.at(0, 1), 2);
    EXPECT_EQ(copy.at(0, 1), 5);
    EXPECT_EQ(copy2.at(0, 1), 6);
}

TEST(Matrix33, At) {
    mat m{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
    EXPECT_EQ(m.at(0, 0), m.data()[0]);
    EXPECT_EQ(m.at(0, 1), m.data()[1]);
    EXPECT_EQ(m.at(0, 2), m.data()[2]);
    EXPECT_EQ(m.at(1, 0), m.data()[3]);
    EXPECT_EQ(m.at(1, 1), m.data()[4]);
    EXPECT_EQ(m.at(1, 2), m.data()[5]);
    EXPECT_EQ(m.at(2, 0), m.data()[6]);
    EXPECT_EQ(m.at(2, 1), m.data()[7]);
    EXPECT_EQ(m.at(2, 2), m.data()[8]);
}

TEST(Matrix33, Plus) {
    mat m2{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
    mat m{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
    m += m2;
    auto m3 = m2 + m2;
    EXPECT_EQ(m.at(0, 0), 2) << "tried accessing (0, 0) = " << m.at(0, 0) << " for matrix " << m;
    EXPECT_EQ(m.at(0, 1), 4) << "tried accessing (0, 1) = " << m.at(0, 1) << " for matrix " << m;
    EXPECT_EQ(m.at(0, 2), 6) << "tried accessing (0, 2) = " << m.at(0, 2) << " for matrix " << m;
    EXPECT_EQ(m.at(1, 0), 8) << "tried accessing (1, 0) = " << m.at(1, 0) << " for matrix " << m;
    EXPECT_EQ(m.at(1, 1), 10) << "tried accessing (1, 1) = " << m.at(1, 1) << " for matrix " << m;
    EXPECT_EQ(m.at(1, 2), 12) << "tried accessing (1, 2) = " << m.at(1, 2) << " for matrix " << m;
    EXPECT_EQ(m.at(2, 0), 14) << "tried accessing (2, 0) = " << m.at(2, 0) << " for matrix " << m;
    EXPECT_EQ(m.at(2, 1), 16) << "tried accessing (2, 1) = " << m.at(2, 1) << " for matrix " << m;
    EXPECT_EQ(m.at(2, 2), 18) << "tried accessing (2, 2) = " << m.at(2, 2) << " for matrix " << m;
    EXPECT_EQ(m3.at(0, 0), 2) << "tried accessing (0, 0) = " << m3.at(0, 0) << " for matrix " << m3;
    EXPECT_EQ(m3.at(0, 1), 4) << "tried accessing (0, 1) = " << m3.at(0, 1) << " for matrix " << m3;
    EXPECT_EQ(m3.at(0, 2), 6) << "tried accessing (0, 2) = " << m3.at(0, 2) << " for matrix " << m3;
    EXPECT_EQ(m3.at(1, 0), 8) << "tried accessing (1, 0) = " << m3.at(1, 0) << " for matrix " << m3;
    EXPECT_EQ(m3.at(1, 1), 10) << "tried accessing (1, 1) = " << m3.at(1, 1) << " for matrix " << m3;
    EXPECT_EQ(m3.at(1, 2), 12) << "tried accessing (1, 2) = " << m3.at(1, 2) << " for matrix " << m3;
    EXPECT_EQ(m3.at(2, 0), 14) << "tried accessing (2, 0) = " << m3.at(2, 0) << " for matrix " << m3;
    EXPECT_EQ(m3.at(2, 1), 16) << "tried accessing (2, 1) = " << m3.at(2, 1) << " for matrix " << m3;
    EXPECT_EQ(m3.at(2, 2), 18) << "tried accessing (2, 2) = " << m3.at(2, 2) << " for matrix " << m3;
    EXPECT_EQ(m, m3);
}

TEST(Matrix33, Scale) {
    mat m{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
    auto m2 = m * 2;
    auto m3 = 2 * m;
    m *= 2;
    EXPECT_EQ(m.at(0, 0), 2);
    EXPECT_EQ(m.at(0, 1), 4);
    EXPECT_EQ(m.at(0, 2), 6);
    EXPECT_EQ(m.at(1, 0), 8);
    EXPECT_EQ(m.at(1, 1), 10);
    EXPECT_EQ(m.at(1, 2), 12);
    EXPECT_EQ(m.at(2, 0), 14);
    EXPECT_EQ(m.at(2, 1), 16);
    EXPECT_EQ(m.at(2, 2), 18);
    EXPECT_EQ(m, m2);
    EXPECT_EQ(m, m3);
}

TEST(Matrix33, EQ) {
    mat m{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
    mat m2{{{1, 2, 3, 4, 5, 6, 7, 8, 9}}};
    mat m3{{{2, 2, 3, 4, 5, 6, 7, 8, 9}}};

    EXPECT_TRUE(m == m);
    EXPECT_TRUE(m == m2);
    EXPECT_FALSE(m == m3);
}

}
