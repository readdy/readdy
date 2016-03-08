//
// Created by clonker on 07.03.16.
//

#include <boost/log/core/core.hpp>
#include "gtest/gtest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    boost::log::core::get()->remove_all_sinks();
    return result;
}

