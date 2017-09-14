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
 *
 *
 * @file KernelConfiguration.h
 * @brief 
 * @author clonker
 * @date 9/13/17
 */
#pragma once

#include <string>
#include <json.hpp>
#include <readdy/common/thread/Config.h>

#include "readdy/common/macros.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(conf)
using json = nlohmann::json;

NAMESPACE_BEGIN(cpu)
struct NeighborList {
    std::string type {"CellDecomposition"};
    std::uint8_t cll_radius {1};
};
void to_json(json &j, const NeighborList &nl);
void from_json(const json &j, NeighborList &nl);

struct ThreadConfig {
    int nThreads {-1};
    util::thread::ThreadMode threadMode {util::thread::ThreadMode::pool};
};
void to_json(json &j, const ThreadConfig &nl);
void from_json(const json &j, ThreadConfig &nl);

struct Configuration {
    NeighborList neighborList;
    ThreadConfig threadConfig;
};
void to_json(json &j, const Configuration &conf);
void from_json(const json &j, Configuration &conf);

NAMESPACE_END(cpu)

NAMESPACE_END(conf)
NAMESPACE_END(readdy)