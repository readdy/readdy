/**
 * << detailed description >>
 *
 * @file TestReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 22.06.16
 */

#include <gtest/gtest.h>
#include <readdy/kernel/singlecpu/reactions/SingleCPUReactions.h>
#include <readdy/model/Kernel.h>
#include <boost/algorithm/string.hpp>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/programs/Programs.h>

namespace {

    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Death;
    using particle_t = readdy::model::Particle;

    class SingleCPUTestReactions : public ::testing::Test {
    protected:
        SingleCPUTestReactions() {
            // if we're in conda
            const char *env = std::getenv("PREFIX");
            std::string pluginDir = "lib/readdy_plugins";
            if (env) {
                auto _env = std::string(env);
                if (!boost::algorithm::ends_with(env, "/")) {
                    _env = _env.append("/");
                }
                pluginDir = _env.append(pluginDir);
            }
            readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory(pluginDir);
        }
    };

    TEST_F(SingleCPUTestReactions, CheckInOutTypesAndPositions) {
        auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
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
            auto fusion = kernel->getReactionFactory().createReaction<fusion_t>("A+B->C", 0, 1, 2, 1, eductDistance, weight1, weight2);
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
            auto fission = kernel->getReactionFactory().createReaction<fission_t>("C->A+B", 2, 0, 1, productDistance, 1, weight1, weight2);
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

    TEST_F(SingleCPUTestReactions, TestDecay) {
        auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
        kernel->getKernelContext().setBoxSize(10, 10, 10);
        kernel->getKernelContext().setTimeStep(1);
        kernel->getKernelContext().setDiffusionConstant("X", .55);
        kernel->getKernelContext().registerDeathReaction("X decay", "X", .05);
        kernel->getKernelContext().registerFissionReaction("X fission", "X", "X", "X", .5, .00);

        auto &&diffuseProgram = kernel->createProgram<readdy::model::programs::DiffuseProgram>();
        auto &&updateModelProgram = kernel->createProgram<readdy::model::programs::UpdateStateModelProgram>();
        auto &&reactionsProgram = kernel->createProgram<readdy::model::programs::DefaultReactionProgram>();

        auto pp_obs = kernel->createObservable<readdy::model::ParticlePositionObservable>();
        pp_obs->setStride(1);
        auto connection = kernel->registerObservable(pp_obs.get());

        const int n_particles = 2000;
        const unsigned int typeId = kernel->getKernelContext().getParticleTypeID("X");
        std::vector<readdy::model::Particle> particlesToBeginWith {n_particles, {0,0,0,typeId}};
        BOOST_LOG_TRIVIAL(debug) << "n_particles="<<particlesToBeginWith.size();
        kernel->getKernelStateModel().addParticles(particlesToBeginWith);

        std::unique_ptr<readdy::model::RandomProvider> rand = std::make_unique<readdy::model::RandomProvider>();

        for(size_t t = 0; t < 1000; t++) {

            diffuseProgram->execute();
            updateModelProgram->configure(t, false);
            updateModelProgram->execute();

            reactionsProgram->execute();
            updateModelProgram->configure(t, false);
            updateModelProgram->execute();

            pp_obs->evaluate();
            BOOST_LOG_TRIVIAL(debug) << "\tcurrently n particles: " << pp_obs->getResult()->size();

        }

        connection.disconnect();
    }

}