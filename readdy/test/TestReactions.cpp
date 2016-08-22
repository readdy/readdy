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
        kernel->registerReaction<readdy::model::reactions::Conversion>("A to B", "A", "B", 0.55);

        {
            // sanity check of operator<< for reactions
            const auto r = kernel->getReactionFactory().createReaction<readdy::model::reactions::Decay>("decay", 0, .1);
            BOOST_LOG_TRIVIAL(debug) << "decay reaction: " << *r;
        }
    }
}