/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


#include <array>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

#include <readdy/common/numeric.h>

class TestTopologyGraphs : public KernelTest {};

/**
 * << detailed description >>
 *
 * @file TestTopologyGraphs.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.03.17
 * @copyright GPL-3
 */

using particle_t = readdy::model::Particle;
using topology_particle_t = readdy::model::TopologyParticle;

using dihedral_bond = readdy::model::top::pot::CosineDihedralPotential;

TEST(TestTopologyGraphs, TestQuadruple) {
    readdy::util::particle_type_quadruple_hasher hasher;
    std::array<readdy::ParticleTypeId, 4> range{1, 2, 3, 4};
    do {
        std::stringstream ss;
        ss << range[0] << ", " << range[1] << ", " << range[2] << ", " << range[3];
        EXPECT_EQ(readdy::util::sortTypeQuadruple(range[0], range[1], range[2], range[3]), std::make_tuple(1, 2, 3, 4))
                            << "failed for range " << ss.str() << ", should always yield a sorted tuple!";
    } while (std::next_permutation(range.begin(), range.end()));
    EXPECT_EQ(hasher(std::make_tuple(1, 2, 3, 4)), hasher(std::make_tuple(4, 3, 2, 1)));
    EXPECT_EQ(hasher(std::make_tuple(1, 3, 2, 4)), hasher(std::make_tuple(4, 2, 3, 1)));
}

TEST(TestTopologyGraphs, TestTriple) {
    readdy::util::particle_type_triple_hasher hasher;
    std::array<readdy::ParticleTypeId, 3> range{1, 2, 3};
    do {
        std::stringstream ss;
        ss << range[0] << ", " << range[1] << ", " << range[2];
        EXPECT_EQ(readdy::util::sortTypeTriple(range[0], range[1], range[2]), std::make_tuple(1, 2, 3))
                            << "failed for range " << ss.str() << ", should always yield a sorted tuple!";
    } while (std::next_permutation(range.begin(), range.end()));
    EXPECT_EQ(hasher(std::make_tuple(1, 2, 3)), hasher(std::make_tuple(3, 2, 1)));
    EXPECT_EQ(hasher(std::make_tuple(2, 1, 3)), hasher(std::make_tuple(3, 1, 2)));
}

TEST(TestTopologyGraphs, TestTuple) {
    readdy::util::particle_type_pair_hasher hasher;
    EXPECT_EQ(hasher(std::make_tuple(1, 2)), hasher(std::make_tuple(2, 1)));
    readdy::util::particle_type_pair_unordered_map<int> map;
    int a = 1, b = 2;
    map[std::make_tuple(a, b)] = 5;
    EXPECT_EQ(map[std::make_tuple(1, 2)], 5);
    EXPECT_EQ(map[std::make_tuple(2, 1)], 5);
}

TEST(TestTopologyGraphs, TestGraphWithIndices) {
    readdy::model::top::graph::Graph graph;
    graph.addVertex(0, 0);
    graph.addVertex(1, 0);
    graph.addEdge(graph.vertices().begin(), ++graph.vertices().begin());
    EXPECT_EQ(graph.vertices().size(), 2);
    EXPECT_EQ(graph.vertexForParticleIndex(0).particleIndex, 0);
    EXPECT_EQ(graph.vertexForParticleIndex(1).particleIndex, 1);
    EXPECT_EQ(graph.vertexForParticleIndex(0).neighbors().size(), 1);
    EXPECT_EQ(graph.vertexForParticleIndex(1).neighbors().size(), 1);
    EXPECT_EQ(graph.vertexForParticleIndex(0).neighbors()[0]->particleIndex, 1);
    EXPECT_EQ(graph.vertexForParticleIndex(1).neighbors()[0]->particleIndex, 0);
    graph.removeParticle(0);
    EXPECT_EQ(graph.vertices().size(), 1);
    EXPECT_EQ(graph.vertexForParticleIndex(1).particleIndex, 1);
    EXPECT_EQ(graph.vertexForParticleIndex(1).neighbors().size(), 0);
}

TEST(TestTopologyGraphs, ConnectedSubComponents) {
    readdy::model::top::graph::Graph graph;
    graph.addVertex(0, 0);
    graph.addVertex(1, 0);
    graph.addVertex(2, 0);

    auto vertex_ref_0 = graph.vertices().begin();
    auto vertex_ref_1 = ++graph.vertices().begin();
    auto vertex_ref_2 = ++(++graph.vertices().begin());

    auto it = graph.vertices().begin();
    auto it_adv = ++graph.vertices().begin();
    graph.addEdge(it, it_adv);

    auto subGraphs = graph.connectedComponentsDestructive();
    ASSERT_EQ(subGraphs.size(), 2);
    {
        ASSERT_EQ(subGraphs[0].vertices().size(), 2);
        ASSERT_EQ(subGraphs[0].vertices().begin(), vertex_ref_0);
        ASSERT_EQ(++subGraphs[0].vertices().begin(), vertex_ref_1);
    }
    {
        ASSERT_EQ(subGraphs[1].vertices().size(), 1);
        ASSERT_EQ(subGraphs[1].vertices().begin(), vertex_ref_2);
    }
}

TEST(TestTopologyGraphs, TestTopologyWithGraph) {
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
    auto &ctx = kernel->context();

    ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.particleTypes().add("Topology B", 1.0, readdy::model::particleflavor::TOPOLOGY);

    ctx.topologyRegistry().configureBondPotential("Topology A", "Topology B", {1.0, 1.0});
    ctx.topologyRegistry().configureBondPotential("Topology A", "Topology A", {1.0, 1.0});
    ctx.topologyRegistry().configureAnglePotential("Topology B", "Topology A", "Topology A", {1.0, 1.0});
    ctx.topologyRegistry().configureTorsionPotential("Topology A", "Topology B", "Topology A", "Topology A",
                                                      {1.0, 1.0, 3.0});

    ctx.boxSize() = {{10, 10, 10}};
    topology_particle_t x_i{-1, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_j{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_k{0, 0, 1, ctx.particleTypes().idOf("Topology B")};
    topology_particle_t x_l{1, .1, 1, ctx.particleTypes().idOf("Topology A")};

    auto top = kernel->stateModel().addTopology(0, {x_i, x_j, x_k, x_l});
    EXPECT_EQ(top->graph().vertices().size(), 4);
    auto it = top->graph().vertices().begin();
    auto it2 = ++top->graph().vertices().begin();
    EXPECT_FALSE(top->graph().isConnected());
    top->graph().addEdge(it++, it2++);
    EXPECT_FALSE(top->graph().isConnected());
    top->graph().addEdge(it++, it2++);
    EXPECT_FALSE(top->graph().isConnected());
    top->graph().addEdge(it++, it2++);
    EXPECT_TRUE(top->graph().isConnected());

    top->graph().addEdge(top->graph().firstVertex(), top->graph().lastVertex());
    EXPECT_TRUE(top->graph().isConnected());
    top->configure();
}

TEST(TestTopologyGraphs, TestFindNTuples) {
    using namespace ::testing;
    readdy::model::top::graph::Graph graph;
    graph.addVertex(0, 0);
    graph.addVertex(1, 0);
    graph.addVertex(2, 0);
    graph.addVertex(3, 0);

    graph.addEdge(graph.firstVertex(), std::next(graph.firstVertex()));
    graph.addEdge(std::next(graph.firstVertex()), std::next(graph.firstVertex(), 2));
    graph.addEdge(std::next(graph.firstVertex(), 2), std::next(graph.firstVertex(), 3));
    graph.addEdge(std::next(graph.firstVertex(), 3), graph.firstVertex());

    auto n_tuples = graph.findNTuples();
    const auto& tuples = std::get<0>(n_tuples);
    const auto& triples = std::get<1>(n_tuples);
    const auto& quadruples = std::get<2>(n_tuples);

    EXPECT_EQ(tuples.size(), 4);
    EXPECT_EQ(triples.size(), 4);
    EXPECT_EQ(quadruples.size(), 4);

    auto a = graph.firstVertex();
    auto b = std::next(graph.firstVertex());
    auto c = std::next(graph.firstVertex(), 2);
    auto d = std::next(graph.firstVertex(), 3);

    // tuples
    ASSERT_THAT(tuples, AnyOf(Contains(std::tie(a, b)), Contains(std::tie(b, a))));
    ASSERT_THAT(tuples, AnyOf(Contains(std::tie(b, c)), Contains(std::tie(c, b))));
    ASSERT_THAT(tuples, AnyOf(Contains(std::tie(c, d)), Contains(std::tie(d, c))));
    ASSERT_THAT(tuples, AnyOf(Contains(std::tie(d, a)), Contains(std::tie(a, d))));

    // triples
    ASSERT_THAT(triples, AnyOf(Contains(std::tie(a, b, c)), Contains(std::tie(c, b, a))));
    ASSERT_THAT(triples, AnyOf(Contains(std::tie(b, c, d)), Contains(std::tie(d, c, b))));
    ASSERT_THAT(triples, AnyOf(Contains(std::tie(c, d, a)), Contains(std::tie(a, d, c))));
    ASSERT_THAT(triples, AnyOf(Contains(std::tie(d, a, b)), Contains(std::tie(b, a, d))));

    // quadruples
    ASSERT_THAT(quadruples, AnyOf(Contains(std::tie(d, a, b, c)), Contains(std::tie(c, b, a, d))));
    ASSERT_THAT(quadruples, AnyOf(Contains(std::tie(a, b, c, d)), Contains(std::tie(d, c, b, a))));
    ASSERT_THAT(quadruples, AnyOf(Contains(std::tie(b, c, d, a)), Contains(std::tie(a, d, c, b))));
    ASSERT_THAT(quadruples, AnyOf(Contains(std::tie(c, d, a, b)), Contains(std::tie(b, a, d, c))));
}


TEST(TestTopologyGraphs, TestFindNTuplesInTriangle) {
    using namespace ::testing;
    readdy::model::top::graph::Graph graph;
    graph.addVertex(0, 0);
    graph.addVertex(1, 0);
    graph.addVertex(2, 0);

    graph.addEdge(graph.firstVertex(), std::next(graph.firstVertex()));
    graph.addEdge(std::next(graph.firstVertex()), std::next(graph.firstVertex(), 2));
    graph.addEdge(std::next(graph.firstVertex(), 2), graph.firstVertex());

    auto n_tuples = graph.findNTuples();
    const auto& tuples = std::get<0>(n_tuples);
    const auto& triples = std::get<1>(n_tuples);
    const auto& quadruples = std::get<2>(n_tuples);

    EXPECT_EQ(tuples.size(), 3);
    EXPECT_EQ(triples.size(), 3);
    EXPECT_EQ(quadruples.size(), 0);

    auto a = graph.firstVertex();
    auto b = std::next(graph.firstVertex());
    auto c = std::next(std::next(graph.firstVertex()));

    // tuples
    ASSERT_THAT(tuples, AnyOf(Contains(std::tie(a, b)), Contains(std::tie(b, a))));
    ASSERT_THAT(tuples, AnyOf(Contains(std::tie(b, c)), Contains(std::tie(c, b))));
    ASSERT_THAT(tuples, AnyOf(Contains(std::tie(c, a)), Contains(std::tie(a, c))));

    // triples
    ASSERT_THAT(triples, AnyOf(Contains(std::tie(a, b, c)), Contains(std::tie(c, b, a))));
    ASSERT_THAT(triples, AnyOf(Contains(std::tie(b, c, a)), Contains(std::tie(a, c, b))));
    ASSERT_THAT(triples, AnyOf(Contains(std::tie(c, a, b)), Contains(std::tie(b, a, c))));
}

TEST_P(TestTopologyGraphs, BondedPotential) {
    auto &ctx = kernel->context();
    auto calculateForces = kernel->actions().calculateForces();
    ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.boxSize() = {{10, 10, 10}};
    ctx.topologyRegistry().configureBondPotential("Topology A", "Topology A", {10, 5});
    topology_particle_t x_i{4, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_j{1, 0, 0, ctx.particleTypes().idOf("Topology A")};
    auto top = kernel->stateModel().addTopology(0, {x_i, x_j});
    top->graph().addEdge(top->graph().vertices().begin(), ++top->graph().vertices().begin());
    top->configure();
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    };

    auto conn = kernel->connectObservable(fObs.get());

    calculateForces->perform();
    kernel->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 2);
    readdy::Vec3 f1{40., 0, 0};
    EXPECT_EQ(collectedForces.at(0), f1);
    readdy::Vec3 f2{-40., 0, 0};
    EXPECT_EQ(collectedForces.at(1), f2);
    EXPECT_EQ(kernel->stateModel().energy(), 40);
}


TEST_P(TestTopologyGraphs, MoreComplicatedAnglePotential) {
    auto &ctx = kernel->context();
    auto calculateForces = kernel->actions().calculateForces();
    ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.boxSize() = {{10, 10, 10}};
    ctx.topologyRegistry().configureBondPotential("Topology A", "Topology A", {0., 1.});
    ctx.topologyRegistry().configureAnglePotential("Topology A", "Topology A", "Topology A",
                                                    {1.0, readdy::util::numeric::pi<readdy::scalar>()});
    topology_particle_t x_i{0.1, 0.1, 0.1, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_j{1.0, 0.0, 0.0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_k{1.0, 0.5, -.3, ctx.particleTypes().idOf("Topology A")};
    auto top = kernel->stateModel().addTopology(0, {x_i, x_j, x_k});
    {
        auto it = top->graph().vertices().begin();
        top->graph().addEdge(std::next(it), it);
        top->graph().addEdge(std::next(it), std::next(it, 2));
    }
    top->configure();
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    };

    auto conn = kernel->connectObservable(fObs.get());

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



TEST_P(TestTopologyGraphs, DihedralPotentialSteeperAngle) {
    auto &ctx = kernel->context();
    auto calculateForces = kernel->actions().calculateForces();
    ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.boxSize() = {{10, 10, 10}};
    topology_particle_t x_i{-1, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_j{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_k{0, 0, 1, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_l{1, 3, 1, ctx.particleTypes().idOf("Topology A")};
    auto top = kernel->stateModel().addTopology(0, {x_i, x_j, x_k, x_l});
    auto it = top->graph().vertices().begin();
    auto it2 = ++top->graph().vertices().begin();
    while(it2 != top->graph().vertices().end()) {
        top->graph().addEdge(it, it2);
        std::advance(it, 1);
        std::advance(it2, 1);
    }
    ctx.topologyRegistry().configureBondPotential("Topology A", "Topology A", {0., 1.});
    kernel->context().topologyRegistry().configureTorsionPotential("Topology A", "Topology A", "Topology A",
                                                                    "Topology A",
                                                                    {1.0, 3, readdy::util::numeric::pi<readdy::scalar>()});
    top->configure();
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    };

    auto conn = kernel->connectObservable(fObs.get());

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

TEST(TestTopologyGraphs, TestAppendParticle) {
    using namespace readdy;
    model::Context context;
    context.topologyRegistry().potentialConfiguration().pairPotentials[std::make_tuple(0, 0)].emplace_back();
    context.topologyRegistry().potentialConfiguration().pairPotentials[std::make_tuple(0, 1)].emplace_back();
    model::top::GraphTopology gt {0, {10, 1, 200}, {0, 0, 0}, context, nullptr};
    {
        auto it = gt.graph().vertices().begin();
        auto it2 = ++gt.graph().vertices().begin();
        gt.graph().addEdge(it, it2);
        gt.graph().addEdge(++it, ++it2);
    }
    gt.configure();
    gt.appendParticle(13, 1, 1, 0);
    gt.configure();

    auto it = std::find_if(gt.graph().vertices().begin(), gt.graph().vertices().end(), [](const model::top::graph::Vertex &v) -> bool {
        return v.particleIndex == 3;
    });
    EXPECT_TRUE(it != gt.graph().vertices().end());
    EXPECT_EQ(it->particleType(), 1);
    auto v2 = std::next(gt.graph().vertices().begin());
    EXPECT_TRUE(gt.graph().containsEdge(it, v2));
}

INSTANTIATE_TEST_CASE_P(TestTopologyGraphsCore, TestTopologyGraphs, ::testing::Values("SingleCPU", "CPU"));

