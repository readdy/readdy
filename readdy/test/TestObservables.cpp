/**
 * << detailed description >>
 *
 * @file TestObservables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 02.05.16
 */

#include <boost/algorithm/string.hpp>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/Simulation.h>
#include <readdy/model/programs/Programs.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

namespace m = readdy::model;

namespace {
    class TestObservables : public KernelTest {

    };

    TEST_P(TestObservables, Foo) {
        readdy::model::_internal::ObservableFactory obsf(kernel.get());
        auto x = obsf.create<m::ParticlePositionObservable>(1);
    }

    TEST_P(TestObservables, TestParticlePositions) {
        const unsigned int n_particles = 100;
        const double diffusionConstant = 1;
        kernel->getKernelContext().setDiffusionConstant("type", diffusionConstant);
        const double timeStep = 1.0;
        kernel->getKernelContext().setTimeStep(timeStep);
        const auto particleTypeId = kernel->getKernelContext().getParticleTypeID("type");
        const auto particles = std::vector<m::Particle>(n_particles, m::Particle(0,0,0, particleTypeId));
        kernel->getKernelStateModel().addParticles(particles);
        auto&& obs = kernel->createObservable<m::ParticlePositionObservable>(3);
        auto &&connection = kernel->connectObservable(obs.get());

        auto&& integrator = kernel->createProgram("Eulerian Brownian dynamics integrator");
        auto&& neighborList = kernel->createProgram<readdy::model::programs::UpdateNeighborList>();
        for(readdy::model::time_step_type t = 0; t < 100; t++) {
            integrator->execute();
            neighborList->execute();
            kernel->evaluateObservables(t);
        }

        const auto& result = obs->getResult();
        const auto&& positions = kernel->getKernelStateModel().getParticlePositions();
        auto it_pos = positions.begin();
        int j = 0;
        for(auto it = result.begin(); it != result.end(); it = std::next(it)) {
            EXPECT_EQ(*it, *it_pos);
            it_pos++;
            ++j;
        }
        EXPECT_TRUE(j == 100);
        connection.disconnect();
    }

    TEST_P(TestObservables, TestCombinerObservable) {
        auto&& o1 = kernel->createObservable<m::ParticlePositionObservable>(1);
        auto&& o2 = kernel->createObservable<m::ParticlePositionObservable>(1);
        auto&& o3 = kernel->createObservable<m::TestCombinerObservable>(o1.get(), o2.get());
        auto&& connection = kernel->connectObservable(o3.get());
        auto&& integrator = kernel->createProgram("Eulerian Brownian dynamics integrator");
        kernel->getKernelStateModel().updateNeighborList();
        for(readdy::model::time_step_type t = 0; t < 100; t++) {
            integrator->execute();
            kernel->getKernelStateModel().updateNeighborList();
        }

        const auto& result = o3->getResult();
        for(auto&& p : result) {
            // todo
        }

        connection.disconnect();
    }

    INSTANTIATE_TEST_CASE_P(TestObservables, TestObservables,
                            ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}

