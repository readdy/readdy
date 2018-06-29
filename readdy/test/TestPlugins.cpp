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
