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
 * Tested against a numerical derivative:
 *
 *   import numpy as np
 *   import numdifftools as nd
 *   import matplotlib.pyplot as plt
 *   from scipy.misc import derivative
 *
 *   def partial_derivative(func, var=0, point=[]):
 *      args = point[:]
 *      def wraps(x):
 *          args[var] = x
 *          return func(*args)
 *      return derivative(wraps, point[var], dx=1e-6)
 *   def partial_derivative_nd(func, var, *args):
 *      args2 = list(args)
 *      def wraps(z):
 *          args2[var] = z
 *          return func(*args2).ravel()
 *      J = nd.Jacobian(wraps, order=8)
 *      return J(args[var].ravel())
 *
 * @file TestTopologies.cpp
 * @brief Tests for topology potentials
 * @author clonker
 * @date 08.02.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/common/numeric.h>
#include <readdy/testing/Utils.h>
#include <readdy/testing/KernelTest.h>

using particle_t = readdy::model::Particle;
using topology_particle_t = readdy::model::TopologyParticle;

using harmonic_bond = readdy::model::top::pot::HarmonicBondPotential;
using angle_bond = readdy::model::top::pot::HarmonicAnglePotential;
using dihedral_bond = readdy::model::top::pot::CosineDihedralPotential;

namespace {


struct TestTopologies : KernelTest {

};

/**
 * k_bonded = 10
 * d_0 = 5
 * def bonded_energy(x_i, x_j):
 *      result = k_bonded * ((np.linalg.norm(x_j - x_i) - d_0) ** 2)
 *      return result
 * x_i = np.array([4.0, 0., 0.])
 * x_j = np.array([1.0, 0., 0.])
 * print bonded_energy(x_i, x_j)
 * partial_derivative_nd(bonded_energy, 0||1, x_i, x_j)
 */
TEST_P(TestTopologies, BondedPotential) {
    auto &ctx = kernel->context();
    ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.boxSize() = {{10, 10, 10}};
    topology_particle_t x_i{4, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_j{1, 0, 0, ctx.particleTypes().idOf("Topology A")};
    auto top = kernel->stateModel().addTopology(0,{x_i, x_j});
    {
        harmonic_bond::bond_configurations bonds;
        bonds.emplace_back(0, 1, 10.0, 5.0);
        top->addBondedPotential<harmonic_bond>(bonds);
    }
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    };

    auto conn = kernel->connectObservable(fObs.get());

    auto calculateForces = kernel->actions().calculateForces();
    calculateForces->perform();
    kernel->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 2);
    readdy::Vec3 f1{40., 0, 0};
    EXPECT_EQ(collectedForces.at(0), f1);
    readdy::Vec3 f2{-40., 0, 0};
    EXPECT_EQ(collectedForces.at(1), f2);
    EXPECT_EQ(kernel->stateModel().energy(), 40);
}

/**
 * theta_0 = np.pi
 * k_angle = 1.0
 * def angle_potential(x_i, x_j, x_k):
 *     x_ji = x_i - x_j
 *     x_jk = x_k - x_j
 *     cos_theta = np.dot(x_ji, x_jk) / (np.linalg.norm(x_ji) * np.linalg.norm(x_jk))
 *     return k_angle * ((np.arccos(cos_theta) - theta_0) ** 2)
 * partial_derivative_nd(angle_potential, 0||1||2, x_i, x_j, x_k)
 */
TEST_P(TestTopologies, AnglePotential) {
    auto &ctx = kernel->context();
    ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.boxSize() = {{10, 10, 10}};
    topology_particle_t x_i{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_j{1, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_k{1, 1, 0, ctx.particleTypes().idOf("Topology A")};
    auto top = kernel->stateModel().addTopology(0, {x_i, x_j, x_k});
    {
        std::vector<angle_bond::angle> angles{{0, 1, 2, 1.0, readdy::util::numeric::pi<readdy::scalar>()}};
        top->addAnglePotential<angle_bond>(std::move(angles));
    }
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    };

    auto conn = kernel->connectObservable(fObs.get());

    auto calculateForces = kernel->actions().calculateForces();
    calculateForces->perform();
    kernel->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 3);
    if(kernel->singlePrecision()) {
        EXPECT_FLOAT_EQ(kernel->stateModel().energy(), static_cast<readdy::scalar>(2.4674011002723395));
    } else {
        EXPECT_DOUBLE_EQ(kernel->stateModel().energy(), static_cast<readdy::scalar>(2.4674011002723395));
    }
    readdy::Vec3 force_x_i{0, -3.14159265, 0};
    readdy::Vec3 force_x_j{-3.14159265, 3.14159265, 0};
    readdy::Vec3 force_x_k{3.14159265, 0., 0.};

    EXPECT_VEC3_NEAR(collectedForces[0], force_x_i, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[1], force_x_j, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[2], force_x_k, 1e-6);
}

TEST_P(TestTopologies, MoreComplicatedAnglePotential) {
    auto &ctx = kernel->context();
    ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.boxSize() = {{10, 10, 10}};
    topology_particle_t x_i{0.1, 0.1, 0.1, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_j{1.0, 0.0, 0.0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_k{1.0, 0.5, -.3, ctx.particleTypes().idOf("Topology A")};
    auto top = kernel->stateModel().addTopology(0, {x_i, x_j, x_k});
    {
        std::vector<angle_bond::angle> angles{{0, 1, 2, 1.0, readdy::util::numeric::pi<readdy::scalar>()}};
        top->addAnglePotential<angle_bond>(std::move(angles));
    }
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    };

    auto conn = kernel->connectObservable(fObs.get());

    auto calculateForces = kernel->actions().calculateForces();
    calculateForces->perform();
    kernel->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 3);
    if(kernel->singlePrecision()) {
        EXPECT_FLOAT_EQ(kernel->stateModel().energy(), static_cast<readdy::scalar>(2.5871244540347655));
    } else {
        EXPECT_DOUBLE_EQ(kernel->stateModel().energy(), static_cast<readdy::scalar>(2.5871244540347655));
    }
    readdy::Vec3 force_x_i{-0.13142034, -3.01536661, 1.83258358};
    readdy::Vec3 force_x_j{-5.32252362, 3.44312692, -1.11964973};
    readdy::Vec3 force_x_k{5.45394396, -0.42776031, -0.71293385};

    EXPECT_VEC3_NEAR(collectedForces[0], force_x_i, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[1], force_x_j, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[2], force_x_k, 1e-6);
}

/**
 * phi_0 = np.pi
 * multiplicity = 3
 * k_dihedral = 1.0
 * def dihedral_potential(x_i, x_j, x_k, x_l):
 *     x_ji = x_i - x_j
 *     x_kj = x_j - x_k
 *     x_kl = x_l - x_k
 *     x_jk = x_k - x_j
 *     m = np.cross(x_ji, x_kj)
 *     n = np.cross(x_kl, x_jk)
 *     cos_phi = np.dot(m, n) / (np.linalg.norm(m) * np.linalg.norm(n))
 *     sin_phi = np.dot(np.cross(m, n), x_jk) / (np.linalg.norm(m) * np.linalg.norm(n) * np.linalg.norm(x_jk))
 *     phi = -np.arctan2(sin_phi, cos_phi)
 *     return k_dihedral * (1 + np.cos(multiplicity * phi - phi_0))
 * partial_derivative_nd(dihedral_potential, 0||1||2||3, x_i, x_j, x_k, x_l) * -1.
 */
TEST_P(TestTopologies, DihedralPotential) {
    auto &ctx = kernel->context();
    ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.boxSize() = {{10, 10, 10}};
    topology_particle_t x_i{-1, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_j{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_k{0, 0, 1, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_l{1, .1, 1, ctx.particleTypes().idOf("Topology A")};
    auto top = kernel->stateModel().addTopology(0, {x_i, x_j, x_k, x_l});
    {
        std::vector<dihedral_bond::dihedral_configuration> dihedrals{{0, 1, 2, 3, 1.0, 3, readdy::util::numeric::pi<readdy::scalar>()}};
        top->addTorsionPotential<dihedral_bond>(dihedrals);
    }
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    };

    auto conn = kernel->connectObservable(fObs.get());
    auto calculateForces = kernel->actions().calculateForces();
    calculateForces->perform();
    kernel->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 4);
    if(kernel->singlePrecision()) {
        EXPECT_FLOAT_EQ(kernel->stateModel().energy(), static_cast<readdy::scalar>(0.044370223263673791));
    } else {
        EXPECT_DOUBLE_EQ(kernel->stateModel().energy(), static_cast<readdy::scalar>(0.044370223263673791));
    }
    readdy::Vec3 force_x_i{0., -0.88371125, 0.};
    readdy::Vec3 force_x_j{0., 0.88371125, 0.};
    readdy::Vec3 force_x_k{-0.08749616, 0.87496163, 0.};
    readdy::Vec3 force_x_l{0.08749616, -0.87496163, 0.};
    EXPECT_VEC3_NEAR(collectedForces[0], force_x_i, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[1], force_x_j, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[2], force_x_k, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[3], force_x_l, 1e-6);
}

TEST_P(TestTopologies, DihedralPotentialSteeperAngle) {
    auto &ctx = kernel->context();
    ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.boxSize() = {{10, 10, 10}};
    topology_particle_t x_i{-1, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_j{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_k{0, 0, 1, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_l{1, 3, 1, ctx.particleTypes().idOf("Topology A")};
    auto top = kernel->stateModel().addTopology(0, {x_i, x_j, x_k, x_l});
    {
        std::vector<dihedral_bond::dihedral_configuration> dihedral{{0, 1, 2, 3, 1.0, 3, readdy::util::numeric::pi<readdy::scalar>()}};
        top->addTorsionPotential(std::make_unique<dihedral_bond>(dihedral));
    }
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    };

    auto conn = kernel->connectObservable(fObs.get());
    auto calculateForces = kernel->actions().calculateForces();
    calculateForces->perform();
    kernel->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 4);
    if(kernel->singlePrecision()) {
        EXPECT_FLOAT_EQ(kernel->stateModel().energy(), static_cast<readdy::scalar>(1.8221921916437787));
    } else {
        EXPECT_DOUBLE_EQ(kernel->stateModel().energy(), static_cast<readdy::scalar>(1.8221921916437787));
    }
    readdy::Vec3 force_x_i{0., 1.70762994, 0.};
    readdy::Vec3 force_x_j{0., -1.70762994, 0.};
    readdy::Vec3 force_x_k{0.51228898, -0.17076299, 0.};
    readdy::Vec3 force_x_l{-0.51228898, 0.17076299, 0.};
    EXPECT_VEC3_NEAR(collectedForces[0], force_x_i, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[1], force_x_j, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[2], force_x_k, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[3], force_x_l, 1e-6);
}

INSTANTIATE_TEST_CASE_P(TestTopologiesCore, TestTopologies, ::testing::Values("SingleCPU", "CPU"));

}