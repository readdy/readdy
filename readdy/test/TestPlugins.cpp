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


//
// Created by clonker on 07.03.16.
//

#include <readdy/plugin/KernelProvider.h>
#include <readdy/testing/KernelMock.h>
#include <readdy/common/string.h>
#include <readdy/common/filesystem.h>

namespace plug = readdy::plugin;

namespace {

TEST(Kernel, LoadingNonexistingPlugin) {
    plug::KernelProvider::getInstance().add("foo", [] { return new readdy::testing::KernelMock("foo"); });
    try {
        plug::KernelProvider::getInstance().create("foo2");
        FAIL() << "Expected NoSuchPluginException!";
    } catch (std::invalid_argument const &ex) {
        SUCCEED() << "invalid argument caught.";
    } catch (...) {
        FAIL() << "Expected NoSuchPluginException!";
    }
}

TEST(Kernel, LoadingExistingPlugin) {
    plug::KernelProvider::getInstance().add("bar", [] { return new readdy::testing::KernelMock("bar"); });
    auto kk_ptr = plug::KernelProvider::getInstance().create("bar");
    EXPECT_STREQ("bar", kk_ptr.get()->name().c_str());
}

TEST(KernelProvider, SanityCheckDefaultDirectory) {
    std::string defaultDirectory = plug::KernelProvider::getInstance().getDefaultKernelDirectory();
    readdy::log::debug("default directory is {}", defaultDirectory);
    SUCCEED();
}

TEST(KernelProvider, TestLoadPluginsFromDirectory) {
    // if we're in conda
    const char *env = std::getenv("CONDA_ENV_PATH");
    if (!env) {
        env = std::getenv("PREFIX");
    }
    std::string pluginDir = "readdy/readdy_plugins";
    if (env) {
        auto _env = std::string(env);
        if (!readdy::util::str::has_suffix(_env, "/")) {
            _env = _env.append("/");
        }
        pluginDir = _env.append(pluginDir);
    }
    plug::KernelProvider::getInstance().loadKernelsFromDirectory(pluginDir);
    readdy::log::debug("current path: {}", readdy::util::fs::current_path());
}

TEST(KernelProvider, TestFoo) {
    auto k = plug::KernelProvider::getInstance().create("SingleCPU");
    auto name = k.get()->name();
    readdy::log::debug("foo name: {}", name);
}

}
