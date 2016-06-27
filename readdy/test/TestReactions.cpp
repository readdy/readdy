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

    TEST(TestReactions, TestReactionFactory) {
        auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
        kernel->getKernelContext().setDiffusionConstant("A", 1.0);
        kernel->getKernelContext().setDiffusionConstant("B", 2.0);
        kernel->getKernelContext().registerConversionReaction("A to B", "A", "B", 0.55);
    }
}