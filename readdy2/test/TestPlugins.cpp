//
// Created by clonker on 07.03.16.
//

#include <plugin/Kernel.h>
#include "gtest/gtest.h"

namespace {
    TEST(Kernel, LoadingPlugins) {
        readdy::plugin::KernelProvider::getInstance();
    }
}
