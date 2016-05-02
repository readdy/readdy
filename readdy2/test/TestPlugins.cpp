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
        auto k = std::make_shared<readdy::model::Kernel>("foo");
        plug::KernelProvider::getInstance().add(std::move(k));
        try {
            plug::KernelProvider::getInstance().get("foo2");
            FAIL() << "Expected NoSuchPluginException!";
        } catch (plug::NoSuchPluginException const &ex) {
            SUCCEED() << "NoSuchPluginException caught.";
        } catch (...) {
            FAIL() << "Expected NoSuchPluginException!";
        }
    }

    TEST(Kernel, LoadingExistingPlugin) {
        auto k = std::make_shared<readdy::model::Kernel>("bar");
        plug::KernelProvider::getInstance().add(std::move(k));
        auto kk_ptr = plug::KernelProvider::getInstance().get("bar");
        EXPECT_STREQ("bar", kk_ptr.get()->getName().c_str());
    }

    TEST(KernelProvider, SanityCheckDefaultDirectory) {
        std::string defaultDirectory = plug::KernelProvider::getInstance().getDefaultKernelDirectory();
        BOOST_LOG_TRIVIAL(debug) << "default directory is " << defaultDirectory;
        SUCCEED();
    }

    TEST(KernelProvider, TestLoadPluginsFromDirectory) {
        // if we're in conda
        const char *env = std::getenv("PREFIX");
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
        auto name = plug::KernelProvider::getInstance().get("SingleCPU").get()->getName();
        BOOST_LOG_TRIVIAL(debug) << "foo name: " << name;
    }

    TEST(KernelProvider, TestTestProgram) {
        auto single_cpu_kernel = plug::KernelProvider::getInstance().get("SingleCPU");
        auto test_program = single_cpu_kernel.get()->createProgram("TestProgram");
        test_program.get()->execute();
    }
}
