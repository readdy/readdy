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
 * @file Utils.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#ifndef READDY_TESTING_UTILS_H
#define READDY_TESTING_UTILS_H

#include <string>
#include <readdy/common/Utils.h>
#include <readdy/common/string.h>
#include <readdy/common/logging.h>

// use this to check floating point equality of readdy's Vec3 object
#define EXPECT_VEC3_EQ(u, v) EXPECT_DOUBLE_EQ((u)[0], (v)[0]); EXPECT_DOUBLE_EQ((u)[1], (v)[1]); EXPECT_DOUBLE_EQ((u)[2], (v)[2])
#define EXPECT_FVEC3_EQ(u, v) EXPECT_FLOAT_EQ((u)[0], (v)[0]); EXPECT_FLOAT_EQ((u)[1], (v)[1]); EXPECT_FLOAT_EQ((u)[2], (v)[2])
#define EXPECT_VEC3_NEAR(u, v, abs_error) EXPECT_NEAR((u)[0], (v)[0], abs_error); EXPECT_NEAR((u)[1], (v)[1], abs_error); EXPECT_NEAR((u)[2], (v)[2], abs_error)

namespace readdy {
namespace testing {

inline std::vector<std::string> getKernelsToTest() {
#ifdef READDY_KERNELS_TO_TEST
    std::vector<std::string> kernels = readdy::util::split(std::string(READDY_KERNELS_TO_TEST), ',');
#else
    std::vector<std::string> kernels{"SingleCPU"};
#endif
    return kernels;
}

inline std::string getPluginsDirectory() {
    // test for several environment variables
    const std::string envs[]{"CONDA_ENV_PATH", "CONDA_PREFIX", "PREFIX"};
    const char *env = nullptr;
    for (auto &&key : envs) {
        env = std::getenv(key.c_str());
        if (env != nullptr) {
            log::trace("Using env-variable for plugin dir prefix {}={}", key, env);
            break;
        }
    }
    std::string pluginDir = "readdy/readdy_plugins";
    if (env != nullptr) {
        auto _env = std::string(env);
        if (!util::str::has_suffix(_env, "/")) {
            _env = _env.append("/");
        }
        pluginDir = _env.append(pluginDir);
    } else {
        log::trace("no environment variables found that indicate plugins dir.");
    }
    return pluginDir;
}
}
}
#endif //READDY_TESTING_UTILS_H
