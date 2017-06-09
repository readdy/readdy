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
#include "gtest/gtest.h"

namespace {

struct Element {
    int val;
    bool deactivated;

    void deactivate() {
        deactivated = true;
    }
};

TEST(Foo, Foo) {
    using namespace readdy;
    util::index_persistent_vector<Element> vec;

    vec.push_back({});
    vec.push_back({});
    vec.push_back({});

    ASSERT_EQ(vec.size(), 3);
}

}