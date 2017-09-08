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
 * @file TestTopologyReactionsExternal.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 16.08.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <gtest/gtest.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/api/Simulation.h>
#include <readdy/model/topologies/Utils.h>
#include <readdy/testing/FloatingPoints.h>

class TestTopologyReactionsExternal : public KernelTest {
protected:
    void SetUp() override {
        auto &ctx = kernel->getKernelContext();
        ctx.particle_types().add("Topology A", 1.0, 1.0, readdy::model::particleflavor::TOPOLOGY);
        ctx.particle_types().add("Topology B", 1.0, 1.0, readdy::model::particleflavor::TOPOLOGY);
        ctx.particle_types().add("Topology Invalid Type", 1.0, 1.0, readdy::model::particleflavor::TOPOLOGY);
        ctx.particle_types().add("A", 1.0, 1.0, readdy::model::particleflavor::NORMAL);
        ctx.particle_types().add("B", 1.0, 1.0, readdy::model::particleflavor::NORMAL);

        ctx.topology_registry().configure_bond_potential("Topology A", "Topology A", {10, 10});
        ctx.topology_registry().configure_bond_potential("Topology A", "Topology B", {10, 10});
        ctx.topology_registry().configure_bond_potential("Topology B", "Topology B", {10, 10});

        ctx.boxSize() = {{10, 10, 10}};
    }
};

namespace {

TEST_P(TestTopologyReactionsExternal, TestTopologyEnzymaticReaction) {
    using namespace readdy;
    auto &ctx = kernel->getKernelContext();
    model::TopologyParticle x_0{c_::zero, c_::zero, c_::zero, ctx.particle_types().id_of("Topology A")};
    kernel->getKernelStateModel().addTopology(0, {x_0});
    kernel->getKernelStateModel().addParticle(
            model::Particle(c_::zero, c_::zero, c_::zero, ctx.particle_types().id_of("A"))
    );
    {
        auto r = kernel->createEnzymaticReaction("TopologyEnzymatic", "Topology A", "A", "B", 1.0, 1.0);
        ctx.reactions().add(std::move(r));
    }
    ctx.configure(false);

    auto particles_beforehand = kernel->getKernelStateModel().getParticles();

    {
        std::size_t nNormalFlavor{0};
        for (const auto &p : particles_beforehand) {
            if (ctx.particle_types().info_of(p.getType()).flavor == readdy::model::particleflavor::NORMAL) {
                ++nNormalFlavor;
            }
        }
        ASSERT_EQ(nNormalFlavor, 1);
    }

    auto nl = kernel->getActionFactory().createAction<readdy::model::actions::UpdateNeighborList>();
    nl->perform();
    auto action = kernel->getActionFactory().createAction<readdy::model::actions::reactions::UncontrolledApproximation>(1.0);
    action->perform();

    auto particles = kernel->getKernelStateModel().getParticles();

    ASSERT_EQ(particles.size(), particles_beforehand.size());
    bool found {false};
    std::size_t nNormalFlavor {0};
    for(const auto& p : particles) {
        if(ctx.particle_types().info_of(p.getType()).flavor == readdy::model::particleflavor::NORMAL) {
            ++nNormalFlavor;
            found |= p.getType() == ctx.particle_types().id_of("B");
        }
    }
    ASSERT_EQ(nNormalFlavor, 1);
    ASSERT_TRUE(found) << "The A particle should have been converted to a B particle";
}

TEST_P(TestTopologyReactionsExternal, TestGetTopologyForParticle) {
    // check that getTopologyForParticle does what it is supposed to do
    using namespace readdy;
    auto &ctx = kernel->getKernelContext();
    model::TopologyParticle x_0{c_::zero, c_::zero, c_::zero, ctx.particle_types().id_of("Topology A")};
    auto toplogy = kernel->getKernelStateModel().addTopology(0, {x_0});
    kernel->getKernelStateModel().addParticle(
            model::Particle(c_::zero, c_::zero, c_::zero, ctx.particle_types().id_of("A"))
    );

    for(auto particle : toplogy->getParticles()) {
        auto returned_top = kernel->getKernelStateModel().getTopologyForParticle(particle);
        ASSERT_EQ(toplogy, returned_top);
    }
}

TEST_P(TestTopologyReactionsExternal, TestGetTopologyForParticleDecay) {
    // check that the particles that are contained in (active) topologies point to their respective topologies, also
    // and especially after topology split reactions
    using namespace readdy;
    Simulation sim;
    auto simKernel = sim.setKernel(std::move(kernel));
    sim.setPeriodicBoundary({{true, true, true}});
    sim.setBoxSize(100, 100, 100);
    sim.registerTopologyType("TA");
    auto topAId = sim.registerParticleType("Topology A", 10., 10., model::particleflavor::TOPOLOGY);
    auto aId = sim.registerParticleType("A", 10., 10.);
    sim.configureTopologyBondPotential("Topology A", "Topology A", 10, 1);

    std::vector<model::TopologyParticle> topologyParticles;
    {
        topologyParticles.reserve(90);
        for (int i = 0; i < 90; ++i) {
            topologyParticles.emplace_back(-49. + i, c_::zero, c_::zero, topAId);
        }
    }
    auto toplogy = sim.addTopology("TA", topologyParticles);
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

    sim.registerStructuralTopologyReaction("TA", r);

    sim.runScheme<api::ReaDDyScheme>(true).evaluateTopologyReactions().configureAndRun(35, 1.);

    log::trace("got n topologies: {}", sim.currentTopologies().size());
    for(auto top : sim.currentTopologies()) {
        for(const auto p : top->getParticles()) {
            ASSERT_EQ(simKernel->getKernelStateModel().getTopologyForParticle(p), top);
        }
    }

}

TEST_P(TestTopologyReactionsExternal, AttachParticle) {
    using namespace readdy;
    Simulation sim;
    auto kernel = sim.setKernel(plugin::KernelProvider::getInstance().create(this->kernel->getName()));
    sim.setPeriodicBoundary({{true, true, true}});
    sim.registerTopologyType("TA");
    sim.setBoxSize(15, 15, 15);
    sim.registerParticleType("middle", c_::zero, c_::one, model::particleflavor::TOPOLOGY);
    sim.registerParticleType("end", c_::zero, c_::one, model::particleflavor::TOPOLOGY);
    sim.registerParticleType("A", c_::zero, c_::one);
    sim.configureTopologyBondPotential("middle", "middle", .00000001, 1);
    sim.configureTopologyBondPotential("middle", "end", .00000001, 1);
    sim.configureTopologyBondPotential("end", "end", .00000001, 1);

    auto top = sim.addTopology("TA", {sim.createTopologyParticle("end", {c_::zero-c_::one, c_::zero, c_::zero}),
                                      sim.createTopologyParticle("middle", {c_::zero, c_::zero, c_::zero}),
                                      sim.createTopologyParticle("end", {c_::zero+c_::one, c_::zero, c_::zero})});
    {
        auto it = top->graph().vertices().begin();
        auto it2 = std::next(top->graph().vertices().begin());
        top->graph().addEdge(it, it2);
        ++it; ++it2;
        top->graph().addEdge(it, it2);
    }

    // register attach reaction that transforms (end, A) -> (middle, end)
    sim.registerSpatialTopologyReaction("attach: TA (end) + (A) -> TA (middle--end)",c_::one, c_::one + c_::half);
    sim.addParticle("A", c_::zero - c_::two, c_::zero, c_::zero);
    sim.addParticle("A", c_::zero - c_::three, c_::zero, c_::zero);
    sim.addParticle("A", c_::zero - c_::four, c_::zero, c_::zero);
    sim.addParticle("A", c_::zero + c_::two, c_::zero, c_::zero);
    sim.addParticle("A", c_::zero + c_::three, c_::zero, c_::zero);
    sim.addParticle("A", c_::zero + c_::four, c_::zero, c_::zero);

    sim.runScheme().evaluateTopologyReactions().configureAndRun(6, 1.);

    const auto& type_registry = kernel->getKernelContext().particle_types();

    EXPECT_TRUE(kernel->getKernelContext().topology_registry().is_spatial_reaction_type("A"));
    EXPECT_TRUE(kernel->getKernelContext().topology_registry().is_spatial_reaction_type("end"));
    EXPECT_FALSE(kernel->getKernelContext().topology_registry().is_spatial_reaction_type("middle"));
    EXPECT_EQ(kernel->getKernelContext().calculateMaxCutoff(), c_::one + c_::half);

    EXPECT_EQ(sim.currentTopologies().size(), 1);
    auto chainTop = sim.currentTopologies().at(0);
    EXPECT_EQ(chainTop->getNParticles(), 3 /*original topology particles*/ + 6 /*attached particles*/);

    auto top_particles = kernel->getKernelStateModel().getParticlesForTopology(*chainTop);

    bool foundEndVertex {false};
    // check that graph is indeed linear
    for(std::size_t idx = 0; idx < chainTop->graph().vertices().size() && !foundEndVertex; ++idx) {
        auto prev_neighbor = std::next(chainTop->graph().vertices().begin(), idx);
        const auto& v_end = *prev_neighbor;
        if(v_end.particleType() == type_registry.id_of("end")) {
            foundEndVertex = true;

            EXPECT_EQ(v_end.neighbors().size(), 1);

            EXPECT_EQ(top_particles.at(idx).getType(), type_registry.id_of("end"));

            using flouble = fp::FloatingPoint<scalar>;
            flouble x_end (top_particles.at(idx).getPos().x);
            flouble y_end (top_particles.at(idx).getPos().y);
            flouble z_end (top_particles.at(idx).getPos().z);
            EXPECT_TRUE(x_end.AlmostEquals(flouble{c_::four}) || x_end.AlmostEquals(flouble{-c_::four}))
                                << "the end particle of our topology sausage should be either at x=4 or x=-4";
            EXPECT_TRUE(y_end.AlmostEquals(flouble{c_::zero})) << "no diffusion going on";
            EXPECT_TRUE(z_end.AlmostEquals(flouble{c_::zero})) << "no diffusion going on";

            auto factor = x_end.AlmostEquals(flouble{c_::four}) ? c_::one : -c_::one;

            // walk along topology sausage, check end particles are always at +-4, the other ones are of type middle
            auto next_neighbor = v_end.neighbors().at(0);
            std::size_t i = 0;
            while(next_neighbor->particleType() != type_registry.id_of("end") && i < 20 /*no endless loop*/) {
                auto next_idx = std::distance(chainTop->graph().vertices().begin(), next_neighbor);
                const auto& next_particle = top_particles.at(static_cast<std::size_t>(next_idx));
                auto predicted_pos = factor*c_::four - factor*(i+1)*c_::one;
                auto actual_pos = next_particle.getPos().x;
                EXPECT_TRUE((flouble(actual_pos).AlmostEquals(flouble(predicted_pos))));
                EXPECT_TRUE((flouble(next_particle.getPos().y)).AlmostEquals(flouble(c_::zero)));
                EXPECT_TRUE((flouble(next_particle.getPos().z)).AlmostEquals(flouble(c_::zero)));
                if(next_neighbor->particleType() == type_registry.id_of("middle")) {
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
    auto &context = kernel->getKernelContext();
    context.topology_registry().add_type("1");
    context.topology_registry().add_type("2");
    context.topology_registry().add_type("3");
    context.topology_registry().add_type("4");
    context.topology_registry().add_type("5");
    context.topology_registry().add_type("6");
    context.topology_registry().add_type("T1");
    context.topology_registry().add_type("T2");
    context.topology_registry().add_type("T3");
    context.topology_registry().add_type("T4");
    context.particle_types().add("p1", 1., 1.);
    context.particle_types().add("p2", 1., 1.);
    context.particle_types().add("p3", 1., 1.);
    context.particle_types().add("p4", 1., 1.);
    context.configure();

    readdy::model::top::reactions::STRParser parser (context.topology_registry());

    {
        auto r = parser.parse("topology-topology fusion type1: T1 (p1) + T2 (p2) -> T3 (p3--p4) [self=true]", 10., 11.);
        ASSERT_EQ(r.mode(), readdy::model::top::reactions::STRMode::TT_FUSION_ALLOW_SELF);
        ASSERT_EQ(r.top_type1(), context.topology_registry().id_of("T1"));
        ASSERT_EQ(r.top_type2(), context.topology_registry().id_of("T2"));
        ASSERT_EQ(r.type1(), context.particle_types().id_of("p1"));
        ASSERT_EQ(r.type2(), context.particle_types().id_of("p2"));
        ASSERT_EQ(r.top_type_to1(), context.topology_registry().id_of("T3"));
        ASSERT_EQ(r.top_type_to2(), readdy::topology_type_empty);
        ASSERT_EQ(r.type_to1(), context.particle_types().id_of("p3"));
        ASSERT_EQ(r.type_to2(), context.particle_types().id_of("p4"));
        ASSERT_EQ(r.name(), "topology-topology fusion type1");
        ASSERT_EQ(r.rate(), 10.);
        ASSERT_EQ(r.radius(), 11.);
    }
    {
        auto r = parser.parse("topology-topology fusion type2: T1 (p1) + T2 (p2) -> T3 (p3--p4)", 10., 11.);
        ASSERT_EQ(r.mode(), readdy::model::top::reactions::STRMode::TT_FUSION);
        ASSERT_EQ(r.top_type1(), context.topology_registry().id_of("T1"));
        ASSERT_EQ(r.top_type2(), context.topology_registry().id_of("T2"));
        ASSERT_EQ(r.type1(), context.particle_types().id_of("p1"));
        ASSERT_EQ(r.type2(), context.particle_types().id_of("p2"));
        ASSERT_EQ(r.top_type_to1(), context.topology_registry().id_of("T3"));
        ASSERT_EQ(r.top_type_to2(), readdy::topology_type_empty);
        ASSERT_EQ(r.type_to1(), context.particle_types().id_of("p3"));
        ASSERT_EQ(r.type_to2(), context.particle_types().id_of("p4"));
        ASSERT_EQ(r.name(), "topology-topology fusion type2");
        ASSERT_EQ(r.rate(), 10.);
        ASSERT_EQ(r.radius(), 11.);
    }
    {
        auto r = parser.parse("topology-topology enzymatic type: T1 (p1) + T2 (p2) -> T3 (p3) + T4 (p4)", 10., 11.);
        ASSERT_EQ(r.mode(), readdy::model::top::reactions::STRMode::TT_ENZYMATIC);
        ASSERT_EQ(r.top_type1(), context.topology_registry().id_of("T1"));
        ASSERT_EQ(r.top_type2(), context.topology_registry().id_of("T2"));
        ASSERT_EQ(r.type1(), context.particle_types().id_of("p1"));
        ASSERT_EQ(r.type2(), context.particle_types().id_of("p2"));
        ASSERT_EQ(r.top_type_to1(), context.topology_registry().id_of("T3"));
        ASSERT_EQ(r.top_type_to2(), context.topology_registry().id_of("T4"));
        ASSERT_EQ(r.type_to1(), context.particle_types().id_of("p3"));
        ASSERT_EQ(r.type_to2(), context.particle_types().id_of("p4"));
        ASSERT_EQ(r.name(), "topology-topology enzymatic type");
        ASSERT_EQ(r.rate(), 10.);
        ASSERT_EQ(r.radius(), 11.);
    }
    {
        auto r = parser.parse("topology-particle fusion type: T1 (p1) + (p2) -> T2 (p3--p4)", 10., 11.);
        ASSERT_EQ(r.mode(), readdy::model::top::reactions::STRMode::TP_FUSION);
        ASSERT_EQ(r.top_type1(), context.topology_registry().id_of("T1"));
        ASSERT_EQ(r.top_type2(), readdy::topology_type_empty);
        ASSERT_EQ(r.type1(), context.particle_types().id_of("p1"));
        ASSERT_EQ(r.type2(), context.particle_types().id_of("p2"));
        ASSERT_EQ(r.top_type_to1(), context.topology_registry().id_of("T2"));
        ASSERT_EQ(r.top_type_to2(), readdy::topology_type_empty);
        ASSERT_EQ(r.type_to1(), context.particle_types().id_of("p3"));
        ASSERT_EQ(r.type_to2(), context.particle_types().id_of("p4"));
        ASSERT_EQ(r.name(), "topology-particle fusion type");
        ASSERT_EQ(r.rate(), 10.);
        ASSERT_EQ(r.radius(), 11.);
    }
    {
        auto r = parser.parse("topology-particle enzymatic type: T1 (p1) + (p2) -> T2 (p3) + (p4)", 10., 11.);
        ASSERT_EQ(r.mode(), readdy::model::top::reactions::STRMode::TP_ENZYMATIC);
        ASSERT_EQ(r.top_type1(), context.topology_registry().id_of("T1"));
        ASSERT_EQ(r.top_type2(), readdy::topology_type_empty);
        ASSERT_EQ(r.type1(), context.particle_types().id_of("p1"));
        ASSERT_EQ(r.type2(), context.particle_types().id_of("p2"));
        ASSERT_EQ(r.top_type_to1(), context.topology_registry().id_of("T2"));
        ASSERT_EQ(r.top_type_to2(), readdy::topology_type_empty);
        ASSERT_EQ(r.type_to1(), context.particle_types().id_of("p3"));
        ASSERT_EQ(r.type_to2(), context.particle_types().id_of("p4"));
        ASSERT_EQ(r.name(), "topology-particle enzymatic type");
        ASSERT_EQ(r.rate(), 10.);
        ASSERT_EQ(r.radius(), 11.);
    }
}

TEST_P(TestTopologyReactionsExternal, AttachTopologies) {
    using namespace readdy;
    Simulation sim;
    auto kernel = sim.setKernel(plugin::KernelProvider::getInstance().create(this->kernel->getName()));
    sim.setPeriodicBoundary({{true, true, true}});
    sim.registerTopologyType("TA");
    sim.setBoxSize(15, 15, 15);
    sim.registerParticleType("middle", c_::zero, c_::one, model::particleflavor::TOPOLOGY);
    sim.registerParticleType("end", c_::zero, c_::one, model::particleflavor::TOPOLOGY);
    sim.configureTopologyBondPotential("middle", "middle", .00000001, 1);
    sim.configureTopologyBondPotential("middle", "end", .00000001, 1);
    sim.configureTopologyBondPotential("end", "end", .00000001, 1);


    for(int i = 0; i < 3; ++i) {
        auto top = sim.addTopology("TA", {sim.createTopologyParticle("end", {c_::zero - c_::four + 3*i, c_::zero, c_::zero}),
                                          sim.createTopologyParticle("middle", {c_::zero - c_::three + 3*i, c_::zero, c_::zero}),
                                          sim.createTopologyParticle("end", {c_::zero - c_::two + 3*i, c_::zero, c_::zero})});
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
    sim.registerSpatialTopologyReaction("merge: TA (end) + TA (end) -> TA (middle--middle)",c_::one, c_::one + c_::half);

    EXPECT_EQ(sim.currentTopologies().size(), 3);
    sim.runScheme().evaluateTopologyReactions().configureAndRun(6, 1.);

    const auto& type_registry = kernel->getKernelContext().particle_types();

    EXPECT_TRUE(kernel->getKernelContext().topology_registry().is_spatial_reaction_type("end"));
    EXPECT_FALSE(kernel->getKernelContext().topology_registry().is_spatial_reaction_type("middle"));
    EXPECT_EQ(kernel->getKernelContext().calculateMaxCutoff(), c_::one + c_::half);

    EXPECT_EQ(sim.currentTopologies().size(), 1);
    auto chainTop = sim.currentTopologies().at(0);
    EXPECT_EQ(chainTop->getNParticles(), 3 /*original topology particles*/ + 6 /*attached particles*/);

    auto top_particles = kernel->getKernelStateModel().getParticlesForTopology(*chainTop);

    bool foundEndVertex {false};
    // check that graph is indeed linear
    for(std::size_t idx = 0; idx < chainTop->graph().vertices().size() && !foundEndVertex; ++idx) {
        auto prev_neighbor = std::next(chainTop->graph().vertices().begin(), idx);
        const auto& v_end = *prev_neighbor;
        if(v_end.particleType() == type_registry.id_of("end")) {
            foundEndVertex = true;

            EXPECT_EQ(v_end.neighbors().size(), 1);

            EXPECT_EQ(top_particles.at(idx).getType(), type_registry.id_of("end"));

            using flouble = fp::FloatingPoint<scalar>;
            flouble x_end (top_particles.at(idx).getPos().x);
            flouble y_end (top_particles.at(idx).getPos().y);
            flouble z_end (top_particles.at(idx).getPos().z);
            EXPECT_TRUE(x_end.AlmostEquals(flouble{c_::four}) || x_end.AlmostEquals(flouble{-c_::four}))
                                << "the end particle of our topology sausage should be either at x=4 or x=-4";
            EXPECT_TRUE(y_end.AlmostEquals(flouble{c_::zero})) << "no diffusion going on";
            EXPECT_TRUE(z_end.AlmostEquals(flouble{c_::zero})) << "no diffusion going on";

            auto factor = x_end.AlmostEquals(flouble{c_::four}) ? c_::one : -c_::one;

            // walk along topology sausage, check end particles are always at +-4, the other ones are of type middle
            auto next_neighbor = v_end.neighbors().at(0);
            std::size_t i = 0;
            while(next_neighbor->particleType() != type_registry.id_of("end") && i < 20 /*no endless loop*/) {
                auto next_idx = std::distance(chainTop->graph().vertices().begin(), next_neighbor);
                const auto& next_particle = top_particles.at(static_cast<std::size_t>(next_idx));
                auto predicted_pos = factor*c_::four - factor*(i+1)*c_::one;
                auto actual_pos = next_particle.getPos().x;
                EXPECT_TRUE((flouble(actual_pos).AlmostEquals(flouble(predicted_pos))));
                EXPECT_TRUE((flouble(next_particle.getPos().y)).AlmostEquals(flouble(c_::zero)));
                EXPECT_TRUE((flouble(next_particle.getPos().z)).AlmostEquals(flouble(c_::zero)));
                if(next_neighbor->particleType() == type_registry.id_of("middle")) {
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
