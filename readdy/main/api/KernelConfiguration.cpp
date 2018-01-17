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
 * @file KernelConfiguration.cpp
 * @brief 
 * @author clonker
 * @date 9/13/17
 */

#include <readdy/api/KernelConfiguration.h>

namespace readdy {
namespace conf {

namespace cpu {
void to_json(json &j, const NeighborList &nl) {
    j = json{{"cll_radius", nl.cll_radius}};
}

void from_json(const json &j, NeighborList &nl) {
    nl.cll_radius = j.at("cll_radius").get<std::uint8_t>();
}

void to_json(json &j, const ThreadConfig &nl) {
    j = json{{"n_threads", nl.nThreads}};
}

void from_json(const json &j, ThreadConfig &nl) {
    if (j.find("n_threads") != j.end()) {
        nl.nThreads = j.at("n_threads").get<int>();
    } else {
        nl.nThreads = readdy_default_n_threads();
    }
}

void to_json(json &j, const Configuration &conf) {
    j = json {{"neighbor_list", conf.neighborList},
              {"thread_config", conf.threadConfig}};
}

void from_json(const json &j, Configuration &conf) {
    if (j.find("neighbor_list") != j.end()) {
        conf.neighborList = j.at("neighbor_list").get<NeighborList>();
    } else {
        conf.neighborList = {};
    }
    if (j.find("thread_config") != j.end()) {
        conf.threadConfig = j.at("thread_config").get<ThreadConfig>();
    } else {
        conf.threadConfig = {};
    }
}
}

void to_json(json &j, const Configuration &conf) {
    j = json {{"CPU", conf.cpu}};
}
void from_json(const json &j, Configuration &conf) {
    if(j.find("CPU") != j.end()) {
        conf.cpu = j.at("CPU").get<cpu::Configuration>();
    }
}

}
}