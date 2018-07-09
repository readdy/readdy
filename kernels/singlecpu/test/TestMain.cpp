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
 * @file SingleCPUTestMain.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.06.16
 */
#include <readdy/testing/Utils.h>
#include <readdy/plugin/KernelProvider.h>
#include "gtest/gtest.h"

int perform_tests(int argc, char **argv) {
    readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory(readdy::testing::getPluginsDirectory());
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

int main(int argc, char **argv) {
    readdy::log::console()->set_level(spdlog::level::warn);
    return perform_tests(argc, argv);
}
