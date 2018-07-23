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
