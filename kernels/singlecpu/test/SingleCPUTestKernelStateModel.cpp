/**
 * << detailed description >>
 *
 * @file SingleCPUTestKernelStateModel.cpp
 * @brief Test the methods that manipulate time-dependent simulation data.
 * @author chrisfroe
 * @date 25.08.16
 */

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/kernel/singlecpu/SingleCPUKernelStateModel.h>
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder2.h>

namespace scpu = readdy::kernel::singlecpu;
namespace scpum = scpu::model;
namespace m = readdy::model;

TEST(KernelStateModel, calculateForces) {
    // need context, setup particle types and potential
    // add particles, update neighborlist, calcforces and check!
    // Thus this depends on particleData, context and neighborlist functioning correct
    scpu::SingleCPUKernel kernel;
    m::KernelContext &ctx = kernel.getKernelContext();
    scpu::SingleCPUKernelStateModel &stateModel = kernel.getKernelStateModel();

    {
        // two A particles with radius 1. -> cutoff 2, distance 1.8 -> r-r_0 = 0.2 -> force = 0.2
        ctx.setDiffusionConstant("A", 1.0);
        ctx.setParticleRadius("A", 1.0);
        ctx.setBoxSize(4.,4.,4.);
        ctx.setPeriodicBoundary(false, false, false);
        scpu::potentials::HarmonicRepulsion pot(&kernel);
        pot.setForceConstant(1.);
        boost::uuids::uuid potId = ctx.registerOrder2Potential(&pot, "A", "A");
        ctx.configure();
        auto typeIdA = ctx.getParticleTypeID("A");
        auto twoParticles = std::vector<m::Particle> {m::Particle(0.,0.,0., typeIdA), m::Particle(0.,0.,1.8, typeIdA)};
        stateModel.addParticles(twoParticles);
        stateModel.updateNeighborList();
        stateModel.calculateForces();
        // check results
        auto data = stateModel.getParticleData();
        auto forcesIt = data->cbegin_forces();
        auto vec0 = m::Vec3(0,0,-0.2);
        auto vec1 = m::Vec3(0,0,0.2);
        std::cout << "force[0] " << *forcesIt << " vec0 " << vec0 <<std::endl;
        std::cout << (*forcesIt == vec0) << std::endl;
        EXPECT_TRUE(*forcesIt == vec0);
        ++forcesIt;
        std::cout << "force[0]" << *forcesIt << std::endl;
        EXPECT_TRUE(*forcesIt == vec1);
        ctx.deregisterPotential(potId);
    }
    //

    EXPECT_TRUE(true);
}