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
 * This header file contains the definitions for objects that can be used to further configure kernels. They all
 * offer (de)serialization from and to json.
 *
 * @file KernelConfiguration.h
 * @brief Definitions for kernel configuration objects
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
/**
 * Struct with configuration attributes for the CPU neighbor list implementations
 */
struct NeighborList {
    /**
     * The radius in the box space of cell-linked lists to consider as neighboring cells. A larger value can 
     * drastically increase memory requirements.
     */
    std::uint8_t cll_radius {1};
};
/**
 * Json serialization of NeighborList config struct
 * @param j the json object
 * @param nl the configurational object
 */
void to_json(json &j, const NeighborList &nl);
/**
 * Json deserialization to NeighborList config struct
 * @param j the json object
 * @param nl the configurational object
 */
void from_json(const json &j, NeighborList &nl);

/**
 * Struct with configuration members that are used to parameterize the threading behavoir of the CPU kernel.
 */
struct ThreadConfig {
    /**
     * Number of threads to use. Per default:
     *     * 4 * n_cores in case of a RELEASE build
     *     * n_cores in case of a DEBUG build
     *     * the value of the environment variable READDY_N_CORES, if set (superseeds the other two options)
     */
    int nThreads {-1};

    int getNThreads() const {
        if(nThreads >= 0) {
            return nThreads;
        }
        return readdy_default_n_threads();
    }
};
/**
 * Json serialization of ThreadConfig
 * @param j the json object
 * @param nl the config
 */
void to_json(json &j, const ThreadConfig &nl);
/**
 * Json deserialization to ThreadConfig
 * @param j the json object
 * @param nl the config
 */
void from_json(const json &j, ThreadConfig &nl);

/**
 * Struct that contains configuration information for the CPU kernel.
 */
struct Configuration {
    /**
     * Configuration of the neighbor list
     */
    NeighborList neighborList {};
    /**
     * Configuration of the threading behavior
     */
    ThreadConfig threadConfig {};
};
/**
 * Json serialization of ThreadConfig
 * @param j the json object
 * @param conf the config
 */
void to_json(json &j, const Configuration &conf);
/**
 * Json deserialization of ThreadConfig
 * @param j the json object
 * @param conf the config
 */
void from_json(const json &j, Configuration &conf);
NAMESPACE_END(cpu)

NAMESPACE_BEGIN(scpu)
/**
 * Struct that contains configuration information for the SingleCPU kernel. Currently empty, as there is nothing to
 * be configured.
 */
struct Configuration {
};
NAMESPACE_END(scpu)

/**
 * Struct that contains configuration information for the SingleCPU as well as the CPU kernel.
 */
struct Configuration {
    /**
     * Configuration for the SingleCPU kernel
     */
    scpu::Configuration scpu;
    /**
     * Configuration for the CPU kernel
     */
    cpu::Configuration cpu;
};
/**
 * Json serialization of the Configuration
 * @param j json object
 * @param conf configuration
 */
void to_json(json &j, const Configuration &conf);
/**
 * Json deserialization of the Configuration
 * @param j json object
 * @param conf configuration
 */
void from_json(const json &j, Configuration &conf);

NAMESPACE_END(conf)
NAMESPACE_END(readdy)
