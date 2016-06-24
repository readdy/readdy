/**
 * << detailed description >>
 *
 * @file TestReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.06.16
 */

#include <boost/algorithm/string.hpp>
#include <readdy/plugin/KernelProvider.h>
#include "gtest/gtest.h"

namespace {

    class TestReactions : public ::testing::Test {
    protected:
        TestReactions() {
            // if we're in conda
            const char *env = std::getenv("CONDA_ENV_PATH");
            std::string pluginDir = "lib/readdy_plugins";
            if (env) {
                auto _env = std::string(env);
                if (!boost::algorithm::ends_with(env, "/")) {
                    _env = _env.append("/");
                }
                pluginDir = _env.append(pluginDir);
            }
            readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory(pluginDir);
        }
    };


    TEST_F(TestReactions, TestReactionFactory) {
        auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
        kernel->getKernelContext().setDiffusionConstant("A", 1.0);
        kernel->getKernelContext().setDiffusionConstant("B", 2.0);
        kernel->getKernelContext().registerConversionReaction("A to B", "A", "B", 0.55);
    }
}