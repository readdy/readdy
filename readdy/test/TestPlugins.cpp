//
// Created by clonker on 07.03.16.
//

#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>
#include <boost/algorithm/string/predicate.hpp>
#include "gtest/gtest.h"

namespace plug = readdy::plugin;

namespace {
    TEST(Kernel, LoadingNonexistingPlugin) {
        plug::KernelProvider::getInstance().add("foo", [] {return new readdy::model::Kernel("foo");});
        try {
            plug::KernelProvider::getInstance().create("foo2");
            FAIL() << "Expected NoSuchPluginException!";
        } catch (plug::NoSuchPluginException const &ex) {
            SUCCEED() << "NoSuchPluginException caught.";
        } catch (...) {
            FAIL() << "Expected NoSuchPluginException!";
        }
    }

    TEST(Kernel, LoadingExistingPlugin) {
        plug::KernelProvider::getInstance().add("bar", [] {return new readdy::model::Kernel("bar");});
        auto kk_ptr = plug::KernelProvider::getInstance().create("bar");
        EXPECT_STREQ("bar", kk_ptr.get()->getName().c_str());
    }

    TEST(KernelProvider, SanityCheckDefaultDirectory) {
        std::string defaultDirectory = plug::KernelProvider::getInstance().getDefaultKernelDirectory();
        BOOST_LOG_TRIVIAL(debug) << "default directory is " << defaultDirectory;
        SUCCEED();
    }

    TEST(KernelProvider, TestLoadPluginsFromDirectory) {
        // if we're in conda
        const char *env = std::getenv("CONDA_ENV_PATH");
        if(!env) {
            env = std::getenv("PREFIX");
        }
        std::string pluginDir = "lib/readdy_plugins";
        if (env) {
            auto _env = std::string(env);
            if (!boost::algorithm::ends_with(env, "/")) {
                _env = _env.append("/");
            }
            pluginDir = _env.append(pluginDir);
        }
        plug::KernelProvider::getInstance().loadKernelsFromDirectory(pluginDir);
        BOOST_LOG_TRIVIAL(debug) << "current path: " << boost::filesystem::current_path().string();
    }

    TEST(KernelProvider, TestFoo) {
        auto k = plug::KernelProvider::getInstance().create("SingleCPU");
        auto name = k.get()->getName();
        BOOST_LOG_TRIVIAL(debug) << "foo name: " << name;
    }

    TEST(KernelProvider, TestTestProgram) {
        auto single_cpu_kernel = plug::KernelProvider::getInstance().create("SingleCPU");
        auto test_program = single_cpu_kernel.get()->createProgram("Test");
        test_program.get()->execute();
    }
}
