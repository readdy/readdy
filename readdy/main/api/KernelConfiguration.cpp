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