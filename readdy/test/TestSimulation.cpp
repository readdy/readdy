#include <readdy/Simulation.h>
#include "gtest/gtest.h"

#if BOOST_OS_MACOS
#include <array>
#endif

#include <boost/algorithm/string.hpp>
#include <readdy/model/Kernel.h>

using namespace readdy;

namespace {

    struct MSDAggregator {

        double msd = 0;
        long T = 0;
        std::vector<readdy::model::Vec3> initialPositions;
        std::shared_ptr<double> result = std::make_shared<double>(0);

        void operator()(readdy::model::ParticlePositionObservable::result_t positions) {
            auto it_init = initialPositions.begin();
            auto it_pos = positions.begin();
            while(it_pos != positions.end()) {
                msd += (*it_init - *it_pos)*(*it_init - *it_pos);
                ++it_init;
                ++it_pos;
            }
            ++T;
            *result = msd/T;
        }
    };

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
        Simulation simulation;
        simulation.setKernel("SingleCPU");
        simulation.setBoxSize(10, 10, 10);
        unsigned int n_particles = 103;
        double diffusionConstant = 1;
        simulation.registerParticleType("type", diffusionConstant, .1);
        for (auto _ = 0; _ < n_particles; ++_) {
            simulation.addParticle(0, 0, 0, "type");
        }
        double timestep = 1;
        MSDAggregator aggregator;
        aggregator.initialPositions = simulation.getAllParticlePositions();
        simulation.registerObservable<readdy::model::ParticlePositionObservable>(aggregator, 1);
        simulation.run(100, timestep);
        auto positions = simulation.getAllParticlePositions();
        double msd = 0;
        for(auto&& position : positions) {
            msd += position*position;
        }
        msd /= positions.size();
        BOOST_LOG_TRIVIAL(debug) << "mean squared displacement: " << msd;
        BOOST_LOG_TRIVIAL(debug) << "mean squared displacement2: " << *aggregator.result;
    }

    TEST_F(TestSimulation, TestObservables) {
        Simulation simulation;
        simulation.setKernel("SingleCPU");
        simulation.setBoxSize(10, 10, 10);
        unsigned int n_particles = 103;
        double diffusionConstant = 1;
        simulation.registerParticleType("type", diffusionConstant, .1);
        for (auto _ = 0; _ < n_particles; ++_) {
            simulation.addParticle(0, 0, 0, "type");
        }
        double timestep = 1;

        int n_callbacks = 0;
        simulation.registerObservable<readdy::model::ParticlePositionObservable>([&n_callbacks](const readdy::model::ParticlePositionObservable::result_t &result) -> void {
            ++n_callbacks;
            EXPECT_EQ(103, result.size());
        }, 1);
        simulation.run(100, timestep);
        EXPECT_EQ(100, n_callbacks);
    }
}