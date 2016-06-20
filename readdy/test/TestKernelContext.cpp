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
#include <readdy/common/Utils.h>
#include <readdy/common/make_unique.h>

namespace m = readdy::model;

namespace {
    struct NOOPPotential : public m::potentials::PotentialOrder2 {
        NOOPPotential() : PotentialOrder2("no op") { }

        virtual double calculateEnergy(const readdy::model::Vec3 &x_i, const readdy::model::Vec3 &x_j) override {return 0;}
        virtual void calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &x_i, const readdy::model::Vec3 &x_j) override {}
        virtual void calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy, const readdy::model::Vec3 &x_i, const readdy::model::Vec3 &x_j) override {}

        virtual NOOPPotential *replicate() const override { return new NOOPPotential(*this); }


    };

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

    TEST(KernelContext, PotentialMap) {
        m::KernelContext ctx;
        auto p1 = std::make_unique<NOOPPotential>();
        ctx.registerOrder2Potential(p1.get(), "a", "b");
        ctx.registerOrder2Potential(p1.get(), "b", "a");
        auto&& vector = ctx.getOrder2Potentials("b", "a");
        EXPECT_EQ(vector.size(), 2);
    }

}