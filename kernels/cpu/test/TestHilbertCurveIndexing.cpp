/**
 * << detailed description >>
 *
 * @file TestHilbertCurveIndexing.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.11.16
 */

#include "gtest/gtest.h"

#include <readdy/kernel/cpu/util/hilbert.h>
#include <readdy/plugin/KernelProvider.h>

namespace {

/**
 * Plan:
 *  0. change data structure to vector of vectors
 *  1. Sort boxes
 *  2. assign particles to boxes, group
 *  3. after some time, periodically, reorganize the data incrementally
 *     (depending on (#(particles ungrouped)/#(particles total))?)
 */

TEST(TestCurveIndexing, Cells) {

    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
    kernel->getKernelContext().setBoxSize(5, 5, 5);
    const auto &simBoxSize = kernel->getKernelContext().getBoxSize();

    // 2x2x2 = 8 boxes
    const std::array<unsigned int, 3> nCells{{2, 2, 2}};
    const std::array<double, 3> cellWidth{{1, 1, 1}};




}

}