//
// Created by clonker on 07.03.16.
//


#include <readdy/testing/Utils.h>
#include <readdy/plugin/KernelProvider.h>
#include <spdlog/spdlog.h>

#include "gtest/gtest.h"

int perform_tests(int argc, char **argv) {
    const auto dir = readdy::testing::getPluginsDirectory();
    readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory(dir);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

int main(int argc, char **argv) {
    //spdlog::set_async_mode(256);
    auto console = spdlog::stdout_color_mt("console");
    console->set_level(spdlog::level::debug);
    console->set_pattern("[          ] [%Y-%m-%d %H:%M:%S] [%t] [%l] %v");

    int result = perform_tests(argc, argv);

    spdlog::drop_all();
    return result;
}
