#include <Simulation.h>
#include "gtest/gtest.h"

using namespace ReaDDy;

namespace {
    class TestSimulation : public ::testing::Test {
    protected:
        Simulation simulation;

        virtual void SetUp() { }

        virtual void TearDown() { }

    };

    TEST_F(TestSimulation, Foo) {
		simulation.setKBT(42);
        EXPECT_EQ(42, simulation.getKBT());
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}