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


/**
 * << detailed description >>
 *
 * @file TestTopologyReactionsExternal.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 16.08.17
 * @copyright GPL-3
 */

#include <gtest/gtest.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/api/Simulation.h>
#include <readdy/model/topologies/Utils.h>
#include <readdy/common/FloatingPoints.h>

class TestTopologyReactionsExternal : public ::testing::TestWithParam<std::string> {
public:
    readdy::Simulation simulation;
protected:

    TestTopologyReactionsExternal() : simulation(GetParam()) {}

    void SetUp() override {
        auto &ctx = simulation.context();
        ctx.topologyRegistry().addType("T");
        ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
        ctx.particleTypes().add("Topology B", 1.0, readdy::model::particleflavor::TOPOLOGY);
        ctx.particleTypes().add("Topology Invalid Type", 1.0, readdy::model::particleflavor::TOPOLOGY);
        ctx.particleTypes().add("A", 1.0, readdy::model::particleflavor::NORMAL);
        ctx.particleTypes().add("B", 1.0, readdy::model::particleflavor::NORMAL);

        ctx.topologyRegistry().configureBondPotential("Topology A", "Topology A", {10, 10});
        ctx.topologyRegistry().configureBondPotential("Topology A", "Topology B", {10, 10});
        ctx.topologyRegistry().configureBondPotential("Topology B", "Topology B", {10, 10});

        ctx.boxSize() = {{10, 10, 10}};
    }
};

namespace {

TEST_P(TestTopologyReactionsExternal, TestTopologyEnzymaticReaction) {
    using namespace readdy;
    auto &ctx = simulation.context();
    model::TopologyParticle x_0{c_::zero, c_::zero, c_::zero, ctx.particleTypes().idOf("Topology A")};
    {
        auto tid = ctx.topologyRegistry().addType("MyType");
        simulation.stateModel().addTopology(tid, {x_0});
    }
    simulation.stateModel().addParticle(
            model::Particle(c_::zero, c_::zero, c_::zero, ctx.particleTypes().idOf("A"))
    );
    ctx.reactions().addEnzymatic("TopologyEnzymatic", "Topology A", "A", "B", 1e16, 1.0);

    auto particles_beforehand = simulation.stateModel().getParticles();

    {
        std::size_t nNormalFlavor{0};
        for (const auto &p : particles_beforehand) {
            if (ctx.particleTypes().infoOf(p.type()).flavor == readdy::model::particleflavor::NORMAL) {
                ++nNormalFlavor;
            }
        }
        ASSERT_EQ(nNormalFlavor, 1);
    }

    auto loop = simulation.createLoop(1.);
    loop.useReactionScheduler("UncontrolledApproximation");
    loop.runInitializeNeighborList();
    loop.runUpdateNeighborList();
    loop.runReactions();

    auto particles = simulation.stateModel().getParticles();

    ASSERT_EQ(particles.size(), particles_beforehand.size());
    bool found {false};
    std::size_t nNormalFlavor {0};
    for(const auto& p : particles) {
        if(ctx.particleTypes().infoOf(p.type()).flavor == readdy::model::particleflavor::NORMAL) {
            ++nNormalFlavor;
            found |= p.type() == ctx.particleTypes().idOf("B");
        }
    }
    ASSERT_EQ(nNormalFlavor, 1);
    ASSERT_TRUE(found) << "The A particle should have been converted to a B particle";
}

TEST_P(TestTopologyReactionsExternal, TestGetTopologyForParticle) {
    // check that getTopologyForParticle does what it is supposed to do
    using namespace readdy;
    auto &ctx = simulation.context();
    model::TopologyParticle x_0{c_::zero, c_::zero, c_::zero, ctx.particleTypes().idOf("Topology A")};
    auto toplogy = simulation.addTopology("T", {x_0});
    simulation.addParticle("A", 0, 0, 0);

    for(auto particle : toplogy->getParticles()) {
        auto returned_top = simulation.stateModel().getTopologyForParticle(particle);
        ASSERT_EQ(toplogy, returned_top);
    }
}

TEST_P(TestTopologyReactionsExternal, TestGetTopologyForParticleDecay) {
    // check that the particles that are contained in (active) topologies point to their respective topologies, also
    // and especially after topology split reactions
    using namespace readdy;
    simulation.context().periodicBoundaryConditions() = {{true, true, true}};
    simulation.context().boxSize() = {{100, 100, 100}};
    simulation.context().topologyRegistry().addType("TA");
    //auto topAId = sim.registerParticleType("Topology A", 10., model::particleflavor::TOPOLOGY);
    //auto aId = sim.registerParticleType("A", 10.);
    auto topAId = simulation.context().particleTypes().idOf("Topology A");
    auto aId = simulation.context().particleTypes().idOf("A");
    simulation.context().topologyRegistry().configureBondPotential("Topology A", "Topology A", {10, 1});

    std::vector<model::TopologyParticle> topologyParticles;
    {
        topologyParticles.reserve(90);
        for (int i = 0; i < 90; ++i) {
            topologyParticles.emplace_back(-49. + i, c_::zero, c_::zero, topAId);
        }
    }
    auto toplogy = simulation.addTopology("TA", topologyParticles);
    {
        auto& graph = toplogy->graph();
        auto it1 = graph.vertices().begin();
        auto it2 = std::next(graph.vertices().begin(), 1);
        while(it2 != graph.vertices().end()) {
            graph.addEdge(it1, it2);
            std::advance(it1, 1);
            std::advance(it2, 1);
        }
    }

    model::top::reactions::StructuralTopologyReaction r {[aId](model::top::GraphTopology& top) {
        model::top::reactions::Recipe recipe (top);
        if(top.getNParticles() > 1) {
            auto rnd = model::rnd::uniform_int(0, static_cast<const int>(top.getNParticles() - 2));
            auto it1 = top.graph().vertices().begin();
            auto it2 = std::next(it1, 1);
            if(rnd > 0) {
                std::advance(it1, rnd);
                std::advance(it2, rnd);
            }
            recipe.removeEdge(it1, it2);
        } else {
            recipe.changeParticleType(top.graph().vertices().begin(), aId);
        }
        return recipe;
    }, .7};

    simulation.context().topologyRegistry().addStructuralReaction("TA", r);

    simulation.run(35, 1.);

    log::trace("got n topologies: {}", simulation.currentTopologies().size());
    for(auto top : simulation.currentTopologies()) {
        for(const auto p : top->getParticles()) {
            ASSERT_EQ(simulation.stateModel().getTopologyForParticle(p), top);
        }
    }

}

TEST_P(TestTopologyReactionsExternal, AttachParticle) {
    using namespace readdy;
    Simulation simulation(this->simulation.selectedKernelType());
    simulation.context().periodicBoundaryConditions() = {{true, true, true}};
    simulation.context().topologyRegistry().addType("TA");
    simulation.context().boxSize() = {{15, 15, 15}};
    simulation.context().particleTypes().add("middle", c_::zero, model::particleflavor::TOPOLOGY);
    simulation.context().particleTypes().add("end", c_::zero, model::particleflavor::TOPOLOGY);
    simulation.context().particleTypes().add("A", c_::zero);
    simulation.context().topologyRegistry().configureBondPotential("middle", "middle", {.00000001, 1});
    simulation.context().topologyRegistry().configureBondPotential("middle", "end", {.00000001, 1});
    simulation.context().topologyRegistry().configureBondPotential("end", "end", {.00000001, 1});

    auto top = simulation.addTopology("TA", {simulation.createTopologyParticle("end", {c_::zero-c_::one, c_::zero, c_::zero}),
                                             simulation.createTopologyParticle("middle", {c_::zero, c_::zero, c_::zero}),
                                             simulation.createTopologyParticle("end", {c_::zero+c_::one, c_::zero, c_::zero})});
    {
        auto it = top->graph().vertices().begin();
        auto it2 = std::next(top->graph().vertices().begin());
        top->graph().addEdge(it, it2);
        ++it; ++it2;
        top->graph().addEdge(it, it2);
    }

    // register attach reaction that transforms (end, A) -> (middle, end)
    simulation.context().topologyRegistry().addSpatialReaction("attach: TA (end) + (A) -> TA (middle--end)",
                                                        1e10, c_::one + c_::half);
    simulation.addParticle("A", c_::zero - c_::two, c_::zero, c_::zero);
    simulation.addParticle("A", c_::zero - c_::three, c_::zero, c_::zero);
    simulation.addParticle("A", c_::zero - c_::four, c_::zero, c_::zero);
    simulation.addParticle("A", c_::zero + c_::two, c_::zero, c_::zero);
    simulation.addParticle("A", c_::zero + c_::three, c_::zero, c_::zero);
    simulation.addParticle("A", c_::zero + c_::four, c_::zero, c_::zero);

    simulation.run(6, 1.);

    const auto& type_registry = simulation.context().particleTypes();

    EXPECT_TRUE(simulation.context().topologyRegistry().isSpatialReactionType("A"));
    EXPECT_TRUE(simulation.context().topologyRegistry().isSpatialReactionType("end"));
    EXPECT_FALSE(simulation.context().topologyRegistry().isSpatialReactionType("middle"));
    EXPECT_EQ(simulation.context().calculateMaxCutoff(), c_::one + c_::half);

    EXPECT_EQ(simulation.currentTopologies().size(), 1);
    auto chainTop = simulation.currentTopologies().at(0);
    EXPECT_EQ(chainTop->getNParticles(), 3 /*original topology particles*/ + 6 /*attached particles*/);

    auto top_particles = simulation.stateModel().getParticlesForTopology(*chainTop);

    bool foundEndVertex {false};
    // check that graph is indeed linear
    for(std::size_t idx = 0; idx < chainTop->graph().vertices().size() && !foundEndVertex; ++idx) {
        auto prev_neighbor = std::next(chainTop->graph().vertices().begin(), idx);
        const auto& v_end = *prev_neighbor;
        if(v_end.particleType() == type_registry.idOf("end")) {
            foundEndVertex = true;

            EXPECT_EQ(v_end.neighbors().size(), 1);

            EXPECT_EQ(top_particles.at(idx).type(), type_registry.idOf("end"));

            using flouble = fp::FloatingPoint<scalar>;
            flouble x_end (top_particles.at(idx).pos().x);
            flouble y_end (top_particles.at(idx).pos().y);
            flouble z_end (top_particles.at(idx).pos().z);
            EXPECT_TRUE(x_end.AlmostEquals(flouble{c_::four}) || x_end.AlmostEquals(flouble{-c_::four}))
                                << "the end particle of our topology sausage should be either at x=4 or x=-4";
            EXPECT_TRUE(y_end.AlmostEquals(flouble{c_::zero})) << "no diffusion going on";
            EXPECT_TRUE(z_end.AlmostEquals(flouble{c_::zero})) << "no diffusion going on";

            auto factor = x_end.AlmostEquals(flouble{c_::four}) ? c_::one : -c_::one;

            // walk along topology sausage, check end particles are always at +-4, the other ones are of type middle
            auto next_neighbor = v_end.neighbors().at(0);
            std::size_t i = 0;
            while(next_neighbor->particleType() != type_registry.idOf("end") && i < 20 /*no endless loop*/) {
                auto next_idx = std::distance(chainTop->graph().vertices().begin(), next_neighbor);
                const auto& next_particle = top_particles.at(static_cast<std::size_t>(next_idx));
                auto predicted_pos = factor*c_::four - factor*(i+1)*c_::one;
                auto actual_pos = next_particle.pos().x;
                EXPECT_TRUE((flouble(actual_pos).AlmostEquals(flouble(predicted_pos))));
                EXPECT_TRUE((flouble(next_particle.pos().y)).AlmostEquals(flouble(c_::zero)));
                EXPECT_TRUE((flouble(next_particle.pos().z)).AlmostEquals(flouble(c_::zero)));
                if(next_neighbor->particleType() == type_registry.idOf("middle")) {
                    EXPECT_EQ(next_neighbor->neighbors().size(), 2);
                    if(next_neighbor->neighbors().at(0) == prev_neighbor) {
                        prev_neighbor = next_neighbor;
                        next_neighbor = next_neighbor->neighbors().at(1);
                    } else {
                        prev_neighbor = next_neighbor;
                        next_neighbor = next_neighbor->neighbors().at(0);
                    }

                } else {
                    EXPECT_EQ(next_neighbor->neighbors().size(), 1);
                }

                ++i;
            }

            EXPECT_EQ(i, 3+6-1-1);
        }
    }
    EXPECT_TRUE(foundEndVertex);
}

TEST_P(TestTopologyReactionsExternal, DefinitionParser) {
    auto &context = simulation.context();
    context.topologyRegistry().addType("1");
    context.topologyRegistry().addType("2");
    context.topologyRegistry().addType("3");
    context.topologyRegistry().addType("4");
    context.topologyRegistry().addType("5");
    context.topologyRegistry().addType("6");
    context.topologyRegistry().addType("T1");
    context.topologyRegistry().addType("T2");
    context.topologyRegistry().addType("T3");
    context.topologyRegistry().addType("T4");
    context.particleTypes().add("p1", 1.);
    context.particleTypes().add("p2", 1.);
    context.particleTypes().add("p3", 1.);
    context.particleTypes().add("p4", 1.);

    readdy::model::top::reactions::STRParser parser (context.topologyRegistry());

    {
        auto r = parser.parse("topology-topology fusion type1: T1 (p1) + T2 (p2) -> T3 (p3--p4) [self=true]", 10., 11.);
        ASSERT_EQ(r.mode(), readdy::model::top::reactions::STRMode::TT_FUSION_ALLOW_SELF);
        ASSERT_EQ(r.top_type1(), context.topologyRegistry().idOf("T1"));
        ASSERT_EQ(r.top_type2(), context.topologyRegistry().idOf("T2"));
        ASSERT_EQ(r.type1(), context.particleTypes().idOf("p1"));
        ASSERT_EQ(r.type2(), context.particleTypes().idOf("p2"));
        ASSERT_EQ(r.top_type_to1(), context.topologyRegistry().idOf("T3"));
        ASSERT_EQ(r.top_type_to2(), readdy::EmptyTopologyId);
        ASSERT_EQ(r.type_to1(), context.particleTypes().idOf("p3"));
        ASSERT_EQ(r.type_to2(), context.particleTypes().idOf("p4"));
        ASSERT_EQ(r.name(), "topology-topology fusion type1");
        ASSERT_EQ(r.rate(), 10.);
        ASSERT_EQ(r.radius(), 11.);
    }
    {
        auto r = parser.parse("topology-topology fusion type2: T1 (p1) + T2 (p2) -> T3 (p3--p4)", 10., 11.);
        ASSERT_EQ(r.mode(), readdy::model::top::reactions::STRMode::TT_FUSION);
        ASSERT_EQ(r.top_type1(), context.topologyRegistry().idOf("T1"));
        ASSERT_EQ(r.top_type2(), context.topologyRegistry().idOf("T2"));
        ASSERT_EQ(r.type1(), context.particleTypes().idOf("p1"));
        ASSERT_EQ(r.type2(), context.particleTypes().idOf("p2"));
        ASSERT_EQ(r.top_type_to1(), context.topologyRegistry().idOf("T3"));
        ASSERT_EQ(r.top_type_to2(), readdy::EmptyTopologyId);
        ASSERT_EQ(r.type_to1(), context.particleTypes().idOf("p3"));
        ASSERT_EQ(r.type_to2(), context.particleTypes().idOf("p4"));
        ASSERT_EQ(r.name(), "topology-topology fusion type2");
        ASSERT_EQ(r.rate(), 10.);
        ASSERT_EQ(r.radius(), 11.);
    }
    {
        auto r = parser.parse("topology-topology enzymatic type: T1 (p1) + T2 (p2) -> T3 (p3) + T4 (p4)", 10., 11.);
        ASSERT_EQ(r.mode(), readdy::model::top::reactions::STRMode::TT_ENZYMATIC);
        ASSERT_EQ(r.top_type1(), context.topologyRegistry().idOf("T1"));
        ASSERT_EQ(r.top_type2(), context.topologyRegistry().idOf("T2"));
        ASSERT_EQ(r.type1(), context.particleTypes().idOf("p1"));
        ASSERT_EQ(r.type2(), context.particleTypes().idOf("p2"));
        ASSERT_EQ(r.top_type_to1(), context.topologyRegistry().idOf("T3"));
        ASSERT_EQ(r.top_type_to2(), context.topologyRegistry().idOf("T4"));
        ASSERT_EQ(r.type_to1(), context.particleTypes().idOf("p3"));
        ASSERT_EQ(r.type_to2(), context.particleTypes().idOf("p4"));
        ASSERT_EQ(r.name(), "topology-topology enzymatic type");
        ASSERT_EQ(r.rate(), 10.);
        ASSERT_EQ(r.radius(), 11.);
    }
    {
        auto r = parser.parse("topology-particle fusion type: T1 (p1) + (p2) -> T2 (p3--p4)", 10., 11.);
        ASSERT_EQ(r.mode(), readdy::model::top::reactions::STRMode::TP_FUSION);
        ASSERT_EQ(r.top_type1(), context.topologyRegistry().idOf("T1"));
        ASSERT_EQ(r.top_type2(), readdy::EmptyTopologyId);
        ASSERT_EQ(r.type1(), context.particleTypes().idOf("p1"));
        ASSERT_EQ(r.type2(), context.particleTypes().idOf("p2"));
        ASSERT_EQ(r.top_type_to1(), context.topologyRegistry().idOf("T2"));
        ASSERT_EQ(r.top_type_to2(), readdy::EmptyTopologyId);
        ASSERT_EQ(r.type_to1(), context.particleTypes().idOf("p3"));
        ASSERT_EQ(r.type_to2(), context.particleTypes().idOf("p4"));
        ASSERT_EQ(r.name(), "topology-particle fusion type");
        ASSERT_EQ(r.rate(), 10.);
        ASSERT_EQ(r.radius(), 11.);
    }
    {
        auto r = parser.parse("topology-particle enzymatic type: T1 (p1) + (p2) -> T2 (p3) + (p4)", 10., 11.);
        ASSERT_EQ(r.mode(), readdy::model::top::reactions::STRMode::TP_ENZYMATIC);
        ASSERT_EQ(r.top_type1(), context.topologyRegistry().idOf("T1"));
        ASSERT_EQ(r.top_type2(), readdy::EmptyTopologyId);
        ASSERT_EQ(r.type1(), context.particleTypes().idOf("p1"));
        ASSERT_EQ(r.type2(), context.particleTypes().idOf("p2"));
        ASSERT_EQ(r.top_type_to1(), context.topologyRegistry().idOf("T2"));
        ASSERT_EQ(r.top_type_to2(), readdy::EmptyTopologyId);
        ASSERT_EQ(r.type_to1(), context.particleTypes().idOf("p3"));
        ASSERT_EQ(r.type_to2(), context.particleTypes().idOf("p4"));
        ASSERT_EQ(r.name(), "topology-particle enzymatic type");
        ASSERT_EQ(r.rate(), 10.);
        ASSERT_EQ(r.radius(), 11.);
    }
}

TEST_P(TestTopologyReactionsExternal, AttachTopologies) {
    using namespace readdy;
    simulation.context().periodicBoundaryConditions() = {{true, true, true}};
    simulation.context().topologyRegistry().addType("TA");
    simulation.context().boxSize() = {{15, 15, 15}};
    simulation.context().particleTypes().add("middle", c_::zero, model::particleflavor::TOPOLOGY);
    simulation.context().particleTypes().add("end", c_::zero, model::particleflavor::TOPOLOGY);
    simulation.context().topologyRegistry().configureBondPotential("middle", "middle", {.00000001, 1});
    simulation.context().topologyRegistry().configureBondPotential("middle", "end", {.00000001, 1});
    simulation.context().topologyRegistry().configureBondPotential("end", "end", {.00000001, 1});

    for(int i = 0; i < 3; ++i) {
        auto top = simulation.addTopology("TA", {simulation.createTopologyParticle("end", {c_::zero - c_::four + 3*i, c_::zero, c_::zero}),
                                                 simulation.createTopologyParticle("middle", {c_::zero - c_::three + 3*i, c_::zero, c_::zero}),
                                                 simulation.createTopologyParticle("end", {c_::zero - c_::two + 3*i, c_::zero, c_::zero})});
        {
            auto it = top->graph().vertices().begin();
            auto it2 = std::next(top->graph().vertices().begin());
            top->graph().addEdge(it, it2);
            ++it;
            ++it2;
            top->graph().addEdge(it, it2);
        }
    }

    // register attach reaction that transforms (end, A) -> (middle, end)
    simulation.context().topologyRegistry().addSpatialReaction("merge: TA (end) + TA (end) -> TA (middle--middle)", 1e3, 1.5);

    EXPECT_EQ(simulation.currentTopologies().size(), 3);
    simulation.run(6, 1.);

    const auto& type_registry = simulation.context().particleTypes();

    EXPECT_TRUE(simulation.context().topologyRegistry().isSpatialReactionType("end"));
    EXPECT_FALSE(simulation.context().topologyRegistry().isSpatialReactionType("middle"));
    EXPECT_EQ(simulation.context().calculateMaxCutoff(), c_::one + c_::half);

    EXPECT_EQ(simulation.currentTopologies().size(), 1);
    auto chainTop = simulation.currentTopologies().at(0);
    EXPECT_EQ(chainTop->getNParticles(), 3 /*original topology particles*/ + 6 /*attached particles*/);

    auto top_particles = simulation.stateModel().getParticlesForTopology(*chainTop);

    bool foundEndVertex {false};
    // check that graph is indeed linear
    for(std::size_t idx = 0; idx < chainTop->graph().vertices().size() && !foundEndVertex; ++idx) {
        auto prev_neighbor = std::next(chainTop->graph().vertices().begin(), idx);
        const auto& v_end = *prev_neighbor;
        if(v_end.particleType() == type_registry.idOf("end")) {
            foundEndVertex = true;

            EXPECT_EQ(v_end.neighbors().size(), 1);

            EXPECT_EQ(top_particles.at(idx).type(), type_registry.idOf("end"));

            using flouble = fp::FloatingPoint<scalar>;
            flouble x_end (top_particles.at(idx).pos().x);
            flouble y_end (top_particles.at(idx).pos().y);
            flouble z_end (top_particles.at(idx).pos().z);
            EXPECT_TRUE(x_end.AlmostEquals(flouble{c_::four}) || x_end.AlmostEquals(flouble{-c_::four}))
                                << "the end particle of our topology sausage should be either at x=4 or x=-4";
            EXPECT_TRUE(y_end.AlmostEquals(flouble{c_::zero})) << "no diffusion going on";
            EXPECT_TRUE(z_end.AlmostEquals(flouble{c_::zero})) << "no diffusion going on";

            auto factor = x_end.AlmostEquals(flouble{c_::four}) ? c_::one : -c_::one;

            // walk along topology sausage, check end particles are always at +-4, the other ones are of type middle
            auto next_neighbor = v_end.neighbors().at(0);
            std::size_t i = 0;
            while(next_neighbor->particleType() != type_registry.idOf("end") && i < 20 /*no endless loop*/) {
                auto next_idx = std::distance(chainTop->graph().vertices().begin(), next_neighbor);
                const auto& next_particle = top_particles.at(static_cast<std::size_t>(next_idx));
                auto predicted_pos = factor*c_::four - factor*(i+1)*c_::one;
                auto actual_pos = next_particle.pos().x;
                EXPECT_TRUE((flouble(actual_pos).AlmostEquals(flouble(predicted_pos))));
                EXPECT_TRUE((flouble(next_particle.pos().y)).AlmostEquals(flouble(c_::zero)));
                EXPECT_TRUE((flouble(next_particle.pos().z)).AlmostEquals(flouble(c_::zero)));
                if(next_neighbor->particleType() == type_registry.idOf("middle")) {
                    EXPECT_EQ(next_neighbor->neighbors().size(), 2);
                    if(next_neighbor->neighbors().at(0) == prev_neighbor) {
                        prev_neighbor = next_neighbor;
                        next_neighbor = next_neighbor->neighbors().at(1);
                    } else {
                        prev_neighbor = next_neighbor;
                        next_neighbor = next_neighbor->neighbors().at(0);
                    }

                } else {
                    EXPECT_EQ(next_neighbor->neighbors().size(), 1);
                }

                ++i;
            }

            EXPECT_EQ(i, 3+6-1-1);
        }
    }
    EXPECT_TRUE(foundEndVertex);
}

INSTANTIATE_TEST_CASE_P(Kernels, TestTopologyReactionsExternal,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}
