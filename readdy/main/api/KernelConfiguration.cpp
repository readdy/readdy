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
    j = json{{"type",       nl.type},
             {"cll_radius", nl.cll_radius}};
}

void from_json(const json &j, NeighborList &nl) {
    nl.type = j.at("type").get<std::string>();
    nl.cll_radius = j.at("cll_radius").get<std::uint8_t>();
}

void to_json(json &j, const ThreadConfig &nl) {
    std::string strMode = [&]() -> std::string {
        switch (nl.threadMode) {
            case util::thread::ThreadMode::inactive:
                return "inactive";
            case util::thread::ThreadMode::pool:
                return "pool";
            case util::thread::ThreadMode::std_thread:
                return "std_thread";
            case util::thread::ThreadMode::std_async:
                return "std_async";
        }
    }();
    j = json{{"n_threads", nl.nThreads},
             {"mode",      strMode}};
}

void from_json(const json &j, ThreadConfig &nl) {
    if (j.find("n_threads") != j.end()) nl.nThreads = j.at("n_threads").get<int>();
    std::string mode = j.find("mode") != j.end() ? j.at("mode").get<std::string>() : "pool";
    if (mode == "inactive") {
        nl.threadMode = util::thread::ThreadMode::inactive;
    } else if (mode == "pool") {
        nl.threadMode = util::thread::ThreadMode::pool;
    } else if (mode == "std_thread") {
        nl.threadMode = util::thread::ThreadMode::std_thread;
    } else if (mode == "std_async") {
        nl.threadMode = util::thread::ThreadMode::std_async;
    } else {
        throw std::invalid_argument("Unknown thread mode: \"" + mode
                                    + "\". Possible modes: inactive, pool, std_thread, std_async.");
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