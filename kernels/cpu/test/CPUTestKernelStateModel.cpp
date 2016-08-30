/**
 * << detailed description >>
 *
 * @file CPUTestKernelStateModel.cpp
 * @brief Test the methods that manipulate time-dependent simulation data for the CPU kernel.
 * @author chrisfroe
 * @date 30.08.16
 */

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder2.h>
#include <readdy/testing/Utils.h>

namespace scpu = readdy::kernel::singlecpu;
namespace cpu = readdy::kernel::cpu;
namespace cpum = cpu::model;
namespace m = readdy::model;

namespace {
    TEST(KernelStateModel, calculateForces) {
        // need context, setup particle types and potential
        // add particles, update neighborlist, calcforces and check!
        // Thus this depends on particleData, context and neighborlist functioning correct
        cpu::CPUKernel kernel;
        m::KernelContext &ctx = kernel.getKernelContext();
        cpu::CPUStateModel &stateModel = kernel.getKernelStateModel();

        {
            // two A particles with radius 1. -> cutoff 2, distance 1.8 -> r-r_0 = 0.2 -> force = 0.2
            ctx.setDiffusionConstant("A", 1.0);
            ctx.setParticleRadius("A", 1.0);
            ctx.setBoxSize(4., 4., 4.);
            ctx.setPeriodicBoundary(false, false, false);
            scpu::potentials::HarmonicRepulsion pot(&kernel);
            pot.setForceConstant(1.);
            boost::uuids::uuid potId = ctx.registerOrder2Potential(&pot, "A", "A");
            ctx.configure();
            auto typeIdA = ctx.getParticleTypeID("A");
            auto twoParticles = std::vector<m::Particle> {m::Particle(0., 0., 0., typeIdA), m::Particle(0., 0., 1.8, typeIdA)};
            stateModel.addParticles(twoParticles);
            stateModel.updateNeighborList();
            stateModel.calculateForces();
            stateModel.calculateForces(); // calculating twice should yield the same result. force and energy must not accumulate
            // check results
            auto data = stateModel.getParticleData();
            auto forcesIt = data->cbegin_forces();
            EXPECT_VEC3_EQ(*forcesIt, m::Vec3(0, 0, -0.2));
            ++forcesIt;
            EXPECT_VEC3_EQ(*forcesIt, m::Vec3(0, 0, 0.2));
            EXPECT_DOUBLE_EQ(stateModel.getEnergy(), 0.02);
            // clean up
            ctx.deregisterPotential(potId);
            stateModel.removeAllParticles();
            auto particles = stateModel.getParticles();
            EXPECT_EQ(particles.size(), 0);
            stateModel.clearNeighborList();
        }

        {
            // several particles without potentials -> forces must all be zero
            ctx.setDiffusionConstant("A", 1.0);
            ctx.setDiffusionConstant("B", 1.0);
            ctx.setParticleRadius("A", 1.0);
            ctx.setParticleRadius("B", 2.0);
            ctx.setBoxSize(4., 4., 4.);
            ctx.setPeriodicBoundary(false, false, false);
            ctx.configure();
            auto typeIdA = ctx.getParticleTypeID("A");
            auto typeIdB = ctx.getParticleTypeID("B");
            auto particlesA = std::vector<m::Particle> {
                    m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0.8, 0, typeIdA), m::Particle(0.2, 0, -0.2, typeIdA)
            };
            auto particlesB = std::vector<m::Particle> {
                    m::Particle(0, 0, 0.1, typeIdB), m::Particle(-1.8, 0, 1.8, typeIdB), m::Particle(0.5, 0, 0, typeIdB)
            };
            stateModel.addParticles(particlesA);
            stateModel.addParticles(particlesB);
            stateModel.updateNeighborList();
            stateModel.calculateForces();
            // check results
            auto data = stateModel.getParticleData();
            auto forcesIt = data->cbegin_forces();
            while (forcesIt != data->cend_forces()) {
                EXPECT_VEC3_EQ(*forcesIt, m::Vec3(0, 0, 0));
                ++forcesIt;
            }
            // clean up
            stateModel.removeAllParticles();
            stateModel.clearNeighborList();
        }

        {
            // similar situation as before but now with repulsion between A and B
            ctx.setDiffusionConstant("A", 1.0);
            ctx.setDiffusionConstant("B", 1.0);
            ctx.setParticleRadius("A", 1.0);
            ctx.setParticleRadius("B", 2.0);
            ctx.setBoxSize(10., 10., 10.);
            ctx.setPeriodicBoundary(true,true,false);
            scpu::potentials::HarmonicRepulsion pot(&kernel);
            pot.setForceConstant(1.);
            boost::uuids::uuid potId = ctx.registerOrder2Potential(&pot, "A", "B");
            ctx.configure();
            auto typeIdA = ctx.getParticleTypeID("A");
            auto typeIdB = ctx.getParticleTypeID("B");
            // There are 6 particles. 0-2 are A particles. 3-5 are B particles.
            // The second B particle is a bit further away
            // -> There are forces between all AB pairs except with the second B particle, namely particle 4
            auto particlesA = std::vector<m::Particle> {
                    m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0.8, 0, typeIdA), m::Particle(0.2, 0, -0.2, typeIdA)
            };
            auto particlesB = std::vector<m::Particle> {
                    m::Particle(0, 0, 0.1, typeIdB), m::Particle(-1.8, -4.0, 1.8, typeIdB), m::Particle(0.5, 0, 0, typeIdB)
            };
            stateModel.addParticles(particlesA);
            stateModel.addParticles(particlesB);
            stateModel.updateNeighborList();
            stateModel.calculateForces();
            // handcalculated expectations
            const double energy03 = 4.205;
            const double energy05 = 3.125;
            const double energy13 = 2.4063226755104354;
            const double energy15 = 2.1148056603830194;
            const double energy23 = 3.4833346173608031;
            const double energy25 = 3.4833346173608031;
            const m::Vec3 force03(0, 0, -2.9);
            const m::Vec3 force05(-2.5, 0, 0);
            const m::Vec3 force13(0, 2.1768336301410027, -0.27210420376762534);
            const m::Vec3 force15(-1.08999682000954, 1.743994912015264, 0);
            const m::Vec3 force23(1.4641005886756873, 0, -2.1961508830135306);
            const m::Vec3 force25(-2.1961508830135306, 0, -1.4641005886756873);
            // check results
            auto data = stateModel.getParticleData();
            auto forcesIt = data->cbegin_forces();
            EXPECT_VEC3_EQ(*forcesIt, force03 + force05) << "force on particle 0 = force03 + force05";
            ++forcesIt;
            EXPECT_VEC3_EQ(*forcesIt, force13 + force15) << "force on particle 1 = force13 + force15";
            ++forcesIt;
            EXPECT_VEC3_EQ(*forcesIt, force23 + force25) << "force on particle 2 = force23 + force25";
            ++forcesIt;
            EXPECT_VEC3_EQ(*forcesIt, (-1. * force03) - force13 - force23) << "force on particle 3 = - force03 - force13 - force23";
            ++forcesIt;
            EXPECT_VEC3_EQ(*forcesIt, m::Vec3(0, 0, 0)) << "force on particle 4 = 0";
            ++forcesIt;
            EXPECT_VEC3_EQ(*forcesIt, (-1. * force05) - force15 - force25) << "force on particle 5 = - force05 - force15 - force25";
            const double totalEnergy = energy03 + energy05 + energy13 + energy15 + energy23 + energy25;
            EXPECT_DOUBLE_EQ(stateModel.getEnergy(), totalEnergy);
            // clean up
            ctx.deregisterPotential(potId);
            stateModel.removeAllParticles();
            stateModel.clearNeighborList();
        }
    }
}