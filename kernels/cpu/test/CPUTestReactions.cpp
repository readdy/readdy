/**
 * << detailed description >>
 *
 * @file TestReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 01.09.16
 */

#include <gtest/gtest.h>
#include <readdy/model/Kernel.h>
#include <boost/algorithm/string.hpp>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/kernel/cpu/programs/Reactions.h>

struct fix_n_threads {
    fix_n_threads(readdy::kernel::cpu::CPUKernel *const kernel, unsigned int n)
            : oldValue(static_cast<unsigned int>(kernel->getNThreads())), kernel(kernel) {
        kernel->setNThreads(n);
    }

    ~fix_n_threads() {
        kernel->setNThreads(oldValue);
    }
private:
    const unsigned int oldValue;
    readdy::kernel::cpu::CPUKernel *const kernel;
};

TEST(CPUTestReactions, CheckInOutTypesAndPositions) {
    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Decay;
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
    kernel->getKernelContext().setPeriodicBoundary(false, false, false);
    kernel->getKernelContext().setBoxSize(100, 100, 100);
    const auto diff = kernel->getKernelContext().getShortestDifferenceFun();
    kernel->getKernelContext().setDiffusionConstant("A", .1); // type id 0
    kernel->getKernelContext().setDiffusionConstant("B", .1); // type id 1
    kernel->getKernelContext().setDiffusionConstant("C", .1); // type id 2

    // test conversion
    {
        auto conversion = kernel->getReactionFactory().createReaction<conversion_t>("A->B", 0, 1, 1);
        particle_t p_A{0, 0, 0, 0};
        particle_t p_out{5, 5, 5, 1};
        conversion->perform(p_A, p_A, p_out, p_out);
        EXPECT_EQ(p_out.getType(), conversion->getTypeTo());
        EXPECT_EQ(p_out.getPos(), p_A.getPos());
    }

    // test fusion
    {
        double eductDistance = .4;
        double weight1 = .3, weight2 = .7;
        auto fusion = kernel->getReactionFactory().createReaction<fusion_t>("A+B->C", 0, 1, 2, 1, eductDistance,
                                                                            weight1, weight2);
        particle_t p_out1{50, 50, 50, 70};
        particle_t p_out2{50, 50, 50, 70};
        particle_t p_A{1, 0, 0, 0};
        particle_t p_B{-1, 0, 0, 1};
        fusion->perform(p_A, p_B, p_out1, p_out1);
        fusion->perform(p_B, p_A, p_out2, p_out2);

        EXPECT_EQ(p_out1.getPos(), p_out2.getPos());
        EXPECT_EQ(p_out1.getType(), p_out2.getType());
        EXPECT_EQ(p_out1.getType(), fusion->getTo());

        EXPECT_EQ(readdy::model::Vec3(.4, 0, 0), p_out1.getPos());
    }

    // fission
    {
        double productDistance = .4;
        double weight1 = .3, weight2 = .7;
        auto fission = kernel->getReactionFactory().createReaction<fission_t>("C->A+B", 2, 0, 1, 1, productDistance,
                                                                              weight1, weight2);
        particle_t p_C{0, 0, 0, 2};
        particle_t p_out1{50, 50, 50, 70};
        particle_t p_out2{50, 50, 50, 70};
        fission->perform(p_C, p_C, p_out1, p_out2);

        EXPECT_EQ(p_out1.getType(), fission->getTo1());
        EXPECT_EQ(p_out2.getType(), fission->getTo2());
        auto p_12 = diff(p_out1.getPos(), p_out2.getPos());
        auto p_12_nondirect = p_out2.getPos() - p_out1.getPos();
        EXPECT_EQ(p_12_nondirect, p_12);
        auto distance = sqrt(p_12 * p_12);
        EXPECT_DOUBLE_EQ(productDistance, distance);
    }

    // enzymatic
    {
        auto enzymatic = kernel->getReactionFactory().createReaction<enzymatic_t>("A+C->B+C", 2, 0, 1, 1, .5);
        particle_t p_A{0, 0, 0, 0};
        particle_t p_C{5, 5, 5, 2};
        {
            particle_t p_out1{50, 50, 50, 70};
            particle_t p_out2{50, 50, 50, 70};
            enzymatic->perform(p_A, p_C, p_out1, p_out2);
            if (p_out1.getType() == enzymatic->getCatalyst()) {
                EXPECT_EQ(enzymatic->getCatalyst(), p_out1.getType());
                EXPECT_EQ(enzymatic->getTo(), p_out2.getType());
                EXPECT_EQ(p_C.getPos(), p_out1.getPos());
                EXPECT_EQ(p_A.getPos(), p_out2.getPos());
            } else {
                EXPECT_EQ(enzymatic->getCatalyst(), p_out2.getType());
                EXPECT_EQ(enzymatic->getTo(), p_out1.getType());
                EXPECT_EQ(p_C.getPos(), p_out2.getPos());
                EXPECT_EQ(p_A.getPos(), p_out1.getPos());
            }
        }
        {
            particle_t p_out1{50, 50, 50, 70};
            particle_t p_out2{50, 50, 50, 70};
            enzymatic->perform(p_C, p_A, p_out1, p_out2);
            if (p_out1.getType() == enzymatic->getCatalyst()) {
                EXPECT_EQ(enzymatic->getCatalyst(), p_out1.getType());
                EXPECT_EQ(enzymatic->getTo(), p_out2.getType());
                EXPECT_EQ(p_C.getPos(), p_out1.getPos());
                EXPECT_EQ(p_A.getPos(), p_out2.getPos());
            } else {
                EXPECT_EQ(enzymatic->getCatalyst(), p_out2.getType());
                EXPECT_EQ(enzymatic->getTo(), p_out1.getType());
                EXPECT_EQ(p_C.getPos(), p_out2.getPos());
                EXPECT_EQ(p_A.getPos(), p_out1.getPos());
            }
        }
    }

}

TEST(CPUTestReactions, TestDecay) {
    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Decay;
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
    kernel->getKernelContext().setBoxSize(10, 10, 10);
    kernel->getKernelContext().setTimeStep(1);
    kernel->getKernelContext().setDiffusionConstant("X", .25);
    kernel->registerReaction<death_t>("X decay", "X", 1);
    kernel->registerReaction<fission_t>("X fission", "X", "X", "X", .5, .3);

    auto &&integrator = kernel->createProgram<readdy::model::programs::EulerBDIntegrator>();
    auto &&forces = kernel->createProgram<readdy::model::programs::CalculateForces>();
    auto &&neighborList = kernel->createProgram<readdy::model::programs::UpdateNeighborList>();
    auto &&reactionsProgram = kernel->createProgram<readdy::model::programs::reactions::GillespieParallel>();

    auto pp_obs = kernel->createObservable<readdy::model::ParticlePositionObservable>(1);
    pp_obs->setCallback([](const readdy::model::ParticlePositionObservable::result_t &t) {
        BOOST_LOG_TRIVIAL(trace) << "got n particles=" << t.size();
    });
    auto connection = kernel->connectObservable(pp_obs.get());

    const int n_particles = 200;
    const unsigned int typeId = kernel->getKernelContext().getParticleTypeID("X");
    std::vector<readdy::model::Particle> particlesToBeginWith{n_particles, {0, 0, 0, typeId}};
    kernel->getKernelStateModel().addParticles(particlesToBeginWith);
    kernel->getKernelContext().configure();
    neighborList->execute();
    for (size_t t = 0; t < 20; t++) {

        forces->execute();
        integrator->execute();
        neighborList->execute();
        reactionsProgram->execute();

        kernel->evaluateObservables(t);

    }

    EXPECT_EQ(0, kernel->getKernelStateModel().getParticlePositions().size());

    connection.disconnect();
}

TEST(CPUTestReactions, TestGillespieParallel) {

    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Decay;
    using particle_t = readdy::model::Particle;
    auto kernel = std::make_unique<readdy::kernel::cpu::CPUKernel>();
    kernel->getKernelContext().setBoxSize(10, 10, 12);
    kernel->getKernelContext().setTimeStep(1);
    kernel->getKernelContext().setPeriodicBoundary(true, true, false);

    kernel->getKernelContext().setDiffusionConstant("A", .25);
    double reactionRadius = 1.0;
    kernel->registerReaction<fusion_t>("annihilation", "A", "A", "A", 1.0, reactionRadius);

    auto &&reactionsProgram = kernel->createProgram<readdy::kernel::cpu::programs::reactions::GillespieParallel>();

    const auto typeId = kernel->getKernelContext().getParticleTypeID("A");

    // this particle goes right into the middle, i.e., into the halo region
    kernel->getKernelStateModel().addParticle({0, 0, 0, typeId});
    // these particles go left and right of this particle into the boxes as problematic ones
    kernel->getKernelStateModel().addParticle({0, 0, -.7, typeId});
    kernel->getKernelStateModel().addParticle({0, 0, .7, typeId});
    // these particles are well inside the boxes and should not be considered problematic
    kernel->getKernelStateModel().addParticle({0, 0, -5, typeId});
    kernel->getKernelStateModel().addParticle({0, 0, 5, typeId});

    kernel->getKernelContext().configure();
    // a box width in z direction of 12 should divide into two boxes of 5x5x6 minus the halo region of width 1.0.
    {
        fix_n_threads n_threads {kernel.get(), 2};
        reactionsProgram->execute();
        EXPECT_EQ(1.0, reactionsProgram->getMaxReactionRadius());
        EXPECT_EQ(6.0, reactionsProgram->getBoxWidth());

    }
}
