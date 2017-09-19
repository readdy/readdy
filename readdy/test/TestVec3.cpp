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
#include <readdy/common/ReaDDyVec3.h>

namespace {

using vec_t = readdy::Vec3;

TEST(Vec3, SizeOfVec3) {
    vec_t vec(0,0,0);
    EXPECT_EQ(sizeof(readdy::scalar)*3, sizeof(vec)) << "a vector should have exactly 24 bytes";
}

}