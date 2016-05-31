/**
 * Testfile containing tests for the KernelContext class.
 *
 * @file TestKernelContext.cpp
 * @brief Testfile for KernelContext.
 * @author clonker
 * @date 19.04.16
 */


#include <readdy/model/KernelContext.h>
#include "gtest/gtest.h"

namespace m = readdy::model;

namespace {
    TEST(KernelContext, SetGetKBT) {
        m::KernelContext ctx;
        ctx.setKBT(42);
        EXPECT_EQ(42, ctx.getKBT());
    }

    TEST(KernelContext, PeriodicBoundary) {
        m::KernelContext ctx;
        ctx.setPeriodicBoundary(true, false, true);
        auto boundary = ctx.getPeriodicBoundary();
        EXPECT_TRUE(boundary[0]);
        EXPECT_FALSE(boundary[1]);
        EXPECT_TRUE(boundary[2]);
    }

    TEST(KernelContext, BoxSize) {
        m::KernelContext ctx;
        ctx.setBoxSize(10, 11, 12);
        auto box_size = ctx.getBoxSize();
        EXPECT_EQ(box_size[0], 10);
        EXPECT_EQ(box_size[1], 11);
        EXPECT_EQ(box_size[2], 12);
    }
}