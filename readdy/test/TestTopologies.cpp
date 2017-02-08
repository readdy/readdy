/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          * 
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file TestTopologies.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 08.02.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/common/numeric.h>
#include <readdy/testing/Utils.h>

using particle_t = readdy::model::Particle;
using topology_particle_t = readdy::model::TopologyParticle;

using harmonic_bond = readdy::model::top::HarmonicBondPotential;
using angle_bond = readdy::model::top::HarmonicAnglePotential;
using dihedral_bond = readdy::model::top::CosineDihedralPotential;

namespace {

TEST(TestTopologies, BondedPotential) {
    auto scpu = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    auto &ctx = scpu->getKernelContext();
    ctx.registerParticleType("Topology A", 1.0, 1.0, particle_t::FLAVOR_TOPOLOGY);
    ctx.setBoxSize(10, 10, 10);
    topology_particle_t x_i{4, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_j{1, 0, 0, ctx.getParticleTypeID("Topology A")};
    auto top = scpu->getKernelStateModel().addTopology({x_i, x_j});
    {
        std::vector<harmonic_bond::Bond> bonds{{0, 1, 5.0, 10.0}};
        top->addBondedPotential(std::make_unique<harmonic_bond>(top, bonds));
    }
    auto fObs = scpu->createObservable<readdy::model::observables::Forces>(1);
    std::vector<readdy::model::Vec3> collectedForces;
    fObs->setCallback([&collectedForces](const readdy::model::observables::Forces::result_t &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    });

    auto conn = scpu->connectObservable(fObs.get());

    ctx.configure();
    scpu->getKernelStateModel().calculateForces();
    scpu->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 2);
    readdy::model::Vec3 f1{40., 0, 0};
    EXPECT_EQ(collectedForces.at(0), f1);
    readdy::model::Vec3 f2{-40., 0, 0};
    EXPECT_EQ(collectedForces.at(1), f2);
    EXPECT_EQ(scpu->getKernelStateModel().getEnergy(), 40.);
}

TEST(TestTopologies, AnglePotential) {
    auto scpu = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    auto &ctx = scpu->getKernelContext();
    ctx.registerParticleType("Topology A", 1.0, 1.0, particle_t::FLAVOR_TOPOLOGY);
    ctx.setBoxSize(10, 10, 10);
    topology_particle_t x_i{0, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_j{1, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_k{1, 1, 0, ctx.getParticleTypeID("Topology A")};
    auto top = scpu->getKernelStateModel().addTopology({x_i, x_j, x_k});
    {
        std::vector<angle_bond::Angle> angles{{0, 1, 2, readdy::util::numeric::pi(), 1.0}};
        top->addAnglePotential(std::make_unique<angle_bond>(top, angles));
    }
    auto fObs = scpu->createObservable<readdy::model::observables::Forces>(1);
    std::vector<readdy::model::Vec3> collectedForces;
    fObs->setCallback([&collectedForces](const readdy::model::observables::Forces::result_t &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    });

    auto conn = scpu->connectObservable(fObs.get());

    ctx.configure();
    scpu->getKernelStateModel().calculateForces();
    scpu->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 3);
    EXPECT_DOUBLE_EQ(scpu->getKernelStateModel().getEnergy(), 2.4674011002723395);
    readdy::model::Vec3 force_x_i{0, -3.14159265, 0};
    readdy::model::Vec3 force_x_j{-3.14159265, 3.14159265, 0};
    readdy::model::Vec3 force_x_k{3.14159265, 0., 0.};

    EXPECT_VEC3_NEAR(collectedForces[0], force_x_i, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[1], force_x_j, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[2], force_x_k, 1e-6);
}

TEST(TestTopologies, DihedralPotential) {
    auto scpu = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    auto &ctx = scpu->getKernelContext();
    ctx.registerParticleType("Topology A", 1.0, 1.0, particle_t::FLAVOR_TOPOLOGY);
    ctx.setBoxSize(10, 10, 10);
    topology_particle_t x_i{-1, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_j{0, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_k{0, 0, 1, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_l{1, .1, 1, ctx.getParticleTypeID("Topology A")};
    auto top = scpu->getKernelStateModel().addTopology({x_i, x_j, x_k, x_l});
    {
        std::vector<dihedral_bond::Dihedral> dihedral{{0, 1, 2, 3, 1.0, 3, readdy::util::numeric::pi()}};
        top->addTorsionPotential(std::make_unique<dihedral_bond>(top, dihedral));
    }
    auto fObs = scpu->createObservable<readdy::model::observables::Forces>(1);
    std::vector<readdy::model::Vec3> collectedForces;
    fObs->setCallback([&collectedForces](const readdy::model::observables::Forces::result_t &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    });

    auto conn = scpu->connectObservable(fObs.get());

    ctx.configure();
    scpu->getKernelStateModel().calculateForces();
    scpu->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 4);
    EXPECT_DOUBLE_EQ(scpu->getKernelStateModel().getEnergy(), 0.044370223263673791);
    readdy::model::Vec3 force_x_i{0., -0.88371125, 0.};
    readdy::model::Vec3 force_x_j{0., 0.88371125, 0.};
    readdy::model::Vec3 force_x_k{-0.08749616, 0.87496163, 0.};
    readdy::model::Vec3 force_x_l{0.08749616, -0.87496163, 0.};
    EXPECT_VEC3_NEAR(collectedForces[0], force_x_i, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[1], force_x_j, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[2], force_x_k, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[3], force_x_l, 1e-6);
}

TEST(TestTopologies, DihedralPotentialSteeperAngle) {
    auto scpu = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    auto &ctx = scpu->getKernelContext();
    ctx.registerParticleType("Topology A", 1.0, 1.0, particle_t::FLAVOR_TOPOLOGY);
    ctx.setBoxSize(10, 10, 10);
    topology_particle_t x_i{-1, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_j{0, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_k{0, 0, 1, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_l{1, 3, 1, ctx.getParticleTypeID("Topology A")};
    auto top = scpu->getKernelStateModel().addTopology({x_i, x_j, x_k, x_l});
    {
        std::vector<dihedral_bond::Dihedral> dihedral{{0, 1, 2, 3, 1.0, 3, readdy::util::numeric::pi()}};
        top->addTorsionPotential(std::make_unique<dihedral_bond>(top, dihedral));
    }
    auto fObs = scpu->createObservable<readdy::model::observables::Forces>(1);
    std::vector<readdy::model::Vec3> collectedForces;
    fObs->setCallback([&collectedForces](const readdy::model::observables::Forces::result_t &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    });

    auto conn = scpu->connectObservable(fObs.get());

    ctx.configure();
    scpu->getKernelStateModel().calculateForces();
    scpu->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 4);
    EXPECT_DOUBLE_EQ(scpu->getKernelStateModel().getEnergy(), 1.8221921916437787);
    readdy::model::Vec3 force_x_i{0., 1.70762994, 0.};
    readdy::model::Vec3 force_x_j{0., -1.70762994, 0.};
    readdy::model::Vec3 force_x_k{0.51228898, -0.17076299, 0.};
    readdy::model::Vec3 force_x_l{-0.51228898, 0.17076299, 0.};
    EXPECT_VEC3_NEAR(collectedForces[0], force_x_i, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[1], force_x_j, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[2], force_x_k, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[3], force_x_l, 1e-6);
}

}