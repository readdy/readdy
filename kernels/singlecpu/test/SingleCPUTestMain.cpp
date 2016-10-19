/**
 * << detailed description >>
 *
 * @file SingleCPUTestMain.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.06.16
 */
#include <readdy/testing/Utils.h>
#include <readdy/plugin/KernelProvider.h>
#include "gtest/gtest.h"

int perform_tests(int argc, char **argv) {
    readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory(readdy::testing::getPluginsDirectory());
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

int main(int argc, char **argv) {
    spdlog::set_sync_mode();
    auto console = spdlog::stdout_color_mt("console");
    console->set_level(spdlog::level::debug);
    console->set_pattern("[          ] [%Y-%m-%d %H:%M:%S] [%t] [%l] %v");

    int result = perform_tests(argc, argv);

    spdlog::drop_all();
    return result;
}