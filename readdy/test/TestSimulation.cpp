#include <readdy/Simulation.h>
#include "gtest/gtest.h"

#if BOOST_OS_MACOS
#include <array>
#endif

#include <boost/algorithm/string.hpp>
#include <readdy/model/Kernel.h>

using namespace readdy;

namespace {
    class TestSimulation : public ::testing::Test {
    protected:
        Simulation simulation;

        TestSimulation() {
            simulation.setKernel("SingleCPU");
        }
    };

    TEST_F(TestSimulation, TestKBT) {
        simulation.setKBT(42);
        EXPECT_EQ(42, simulation.getKBT());
    }

    TEST_F(TestSimulation, TestPeriodicBdry) {
        simulation.setPeriodicBoundary({true, false, true});
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

    TEST_F(TestSimulation, TestMeanSquaredDisplacement) {
        simulation.setBoxSize(10, 10, 10);
        unsigned int n_particles = 103;
        double diffusionConstant = 1;
        simulation.registerParticleType("type", diffusionConstant);
        for (auto _ = 0; _ < n_particles; ++_) {
            simulation.addParticle(0, 0, 0, "type");
        }
        double timestep = 1;
        simulation.run(100, timestep);
        auto positions = simulation.getParticlePositions();
        double msd = 0;
        for(auto&& position : positions) {
            msd += position*position;
        }
        msd /= positions.size();
        BOOST_LOG_TRIVIAL(debug) << "mean squared displacement: " << msd;
    }

    TEST_F(TestSimulation, TestObservables) {
        Simulation simulation;
        simulation.setKernel("SingleCPU");
        simulation.setBoxSize(10, 10, 10);
        unsigned int n_particles = 103;
        double diffusionConstant = 1;
        simulation.registerParticleType("type", diffusionConstant);
        for (auto _ = 0; _ < n_particles; ++_) {
            simulation.addParticle(0, 0, 0, "type");
        }
        double timestep = 1;

        int n_callbacks = 0;
        simulation.registerObservable<readdy::model::ParticlePositionObservable>(1, [&n_callbacks](const readdy::model::ParticlePositionObservable::result_t &result) -> void {
            ++n_callbacks;
            EXPECT_EQ(103, result.size());
        });
        simulation.run(100, timestep);
        EXPECT_EQ(100, n_callbacks);
    }
}