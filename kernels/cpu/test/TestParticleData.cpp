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
 * @file TestFoo.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#include <readdy/common/logging.h>
#include <readdy/kernel/cpu/model/CPUParticleData.h>
#include "gtest/gtest.h"

namespace {
    TEST(TestParticleData, TestEntryBytesize) {
        readdy::model::Particle p {0, 0, 0, 0};
        readdy::kernel::cpu::model::CPUParticleData::Entry entry {p};
        EXPECT_EQ(72, sizeof(entry)) << "an entry should have exactly 72 bytes";
    }
}