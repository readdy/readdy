#include <readdy/Simulation.h>
#include "gtest/gtest.h"

using namespace readdy;

namespace {
    class TestSimulation : public ::testing::Test {
    protected:
        Simulation simulation;

        virtual void SetUp() { }

        virtual void TearDown() { }

    };

    TEST_F(TestSimulation, TestKBT) {
		simulation.setKBT(42);
        EXPECT_EQ(42, simulation.getKBT());
    }

    TEST_F(TestSimulation, TestPeriodicBdry) {
        simulation.setPeriodicBoundary(true, false, true);
        auto boundary = simulation.getPeriodicBoundary();
        EXPECT_TRUE(boundary[0]);
        EXPECT_FALSE(boundary[1]);
        EXPECT_TRUE(boundary[2]);
    }

    TEST_F(TestSimulation, TestBoxSize) {
        simulation.setBoxSize(10, 11, 12);
        auto box_size = simulation.getBoxSize();
        EXPECT_EQ(box_size[0], 10);
        EXPECT_EQ(box_size[1], 11);
        EXPECT_EQ(box_size[2], 12);
    }
}