//
// Created by clonker on 07.03.16.
//

#include <readdy/plugin/Kernel.h>
#include "gtest/gtest.h"

namespace plug = readdy::plugin;

namespace {
    TEST(Kernel, LoadingNonexistingPlugin) {
        plug::Kernel k("foo");
        plug::KernelProvider::getInstance().add(k);
        try{
            plug::KernelProvider::getInstance().get("foo2");
            FAIL() << "Expected NoSuchPluginException!";
        } catch(plug::NoSuchPluginException const &ex) {
            SUCCEED() << "NoSuchPluginException caught.";
        } catch(...) {
            FAIL() << "Expected NoSuchPluginException!";
        }
    }

    TEST(Kernel, LoadingExistingPlugin) {
        plug::Kernel k("bar");
        plug::KernelProvider::getInstance().add(k);
        auto kk_ptr = plug::KernelProvider::getInstance().get("bar");
        EXPECT_STREQ("bar", kk_ptr.get()->getName().c_str());
    }

    TEST(KernelProvider, SanityCheckDefaultDirectory) {
        std::string defaultDirectory = plug::KernelProvider::getInstance().getDefaultKernelDirectory();
        BOOST_LOG_TRIVIAL(debug) << "default directory is " << defaultDirectory;
        SUCCEED();
    }

    TEST(KernelProvider, TestLoadPluginsFromDirectory) {
        plug::KernelProvider::getInstance().loadKernelsFromDirectory("out/lib/plugins");
        BOOST_LOG_TRIVIAL(debug) << "foo";
        std::cout << "refcount == " << plug::KernelProvider::getInstance().get("SingleCPU").use_count() << std::endl;
    }

    TEST(KernelProvider, TestFoo) {
        auto name = plug::KernelProvider::getInstance().get("SingleCPU").get()->getName();
        BOOST_LOG_TRIVIAL(debug) << "foo name: " << name;
    }
}
