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
 * @file TestTopologyReactionsExternal.cpp
 * @author clonker
 * @date 16.08.17
 * @copyright BSD-3
 */

#include <catch2/catch.hpp>

#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/api/Simulation.h>
#include <readdy/common/FloatingPoints.h>

using namespace readdy;
using namespace readdytesting::kernel;

TEMPLATE_TEST_CASE("Test topology reactions external", "[topologies]", SingleCPU, CPU) {
    readdy::model::Context ctx;

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

    SECTION("Enzymatic reaction") {
        ctx.reactions().addEnzymatic("TopologyEnzymatic", "Topology A", "A", "B", 1e16, 1.0);

        Simulation simulation {create<TestType>(), ctx};

        model::Particle x_0{0., 0., 0., simulation.context().particleTypes().idOf("Topology A")};
        {
            auto tid = ctx.topologyRegistry().addType("MyType");
            simulation.stateModel().addTopology(tid, {x_0});
        }

        simulation.stateModel().addParticle(
                model::Particle(0., 0., 0., ctx.particleTypes().idOf("A"))
        );
        auto particles_beforehand = simulation.stateModel().getParticles();

        {
            std::size_t nNormalFlavor{0};
            for (const auto &p : particles_beforehand) {
                if (ctx.particleTypes().infoOf(p.type()).flavor == readdy::model::particleflavor::NORMAL) {
                    ++nNormalFlavor;
                }
            }
            REQUIRE(nNormalFlavor == 1);
        }

        auto loop = simulation.createLoop(1.);
        loop.useReactionScheduler("UncontrolledApproximation");
        loop.runInitializeNeighborList();
        loop.runUpdateNeighborList();
        loop.runReactions();

        auto particles = simulation.stateModel().getParticles();

        REQUIRE(particles.size() == particles_beforehand.size());
        bool foundBParticle{false};
        std::size_t nNormalFlavor{0};
        for (const auto &p : particles) {
            if (ctx.particleTypes().infoOf(p.type()).flavor == readdy::model::particleflavor::NORMAL) {
                ++nNormalFlavor;
                foundBParticle |= p.type() == ctx.particleTypes().idOf("B");
            }
        }
        REQUIRE(nNormalFlavor == 1);
        // The A particle should have been converted to a B particle
        REQUIRE(foundBParticle);
    }

    SECTION("Get topology for particle") {
        Simulation simulation {create<TestType>(), ctx};
        model::Particle x_0{0., 0., 0., ctx.particleTypes().idOf("Topology A")};
        auto toplogy = simulation.addTopology("T", {x_0});
        simulation.addParticle("A", 0, 0, 0);

        for(auto particle : toplogy->particleIndices()) {
            auto returned_top = simulation.stateModel().getTopologyForParticle(particle);
            REQUIRE(toplogy == returned_top);
        }
    }

    SECTION("Get topology for particle decay") {
        // check that the particles that are contained in (active) topologies point to their respective topologies, also
        // and especially after topology split reactions
        using namespace readdy;
        
        ctx.periodicBoundaryConditions() = {{true, true, true}};
        ctx.boxSize() = {{100, 100, 100}};
        ctx.topologyRegistry().addType("TA");
        //auto topAId = sim.registerParticleType("Topology A", 10., model::particleflavor::TOPOLOGY);
        //auto aId = sim.registerParticleType("A", 10.);
        auto topAId = ctx.particleTypes().idOf("Topology A");
        auto aId = ctx.particleTypes().idOf("A");
        ctx.topologyRegistry().configureBondPotential("Topology A", "Topology A", {10, 1});

        model::top::reactions::StructuralTopologyReaction r {"r", [aId](model::top::GraphTopology& top) {
            model::top::reactions::Recipe recipe (top);
            if(top.nParticles() > 1) {
                auto rnd = model::rnd::uniform_int(0, static_cast<const int>(top.nParticles() - 2));
                auto it1 = top.graph().begin();
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
        ctx.topologyRegistry().addStructuralReaction("TA", r);

        Simulation simulation {create<TestType>(), ctx};
        
        std::vector<model::TopologyParticle> topologyParticles;
        {
            topologyParticles.reserve(90);
            for (int i = 0; i < 90; ++i) {
                topologyParticles.emplace_back(-49. + i, 0., 0., topAId);
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

        simulation.run(35, 1.);

        log::trace("got n topologies: {}", simulation.currentTopologies().size());
        for(auto top : simulation.currentTopologies()) {
            for(const auto p : top->getParticles()) {
                REQUIRE(simulation.stateModel().getTopologyForParticle(p) == top);
            }
        }
    }

    SECTION("The parser") {
        
        ctx.topologyRegistry().addType("1");
        ctx.topologyRegistry().addType("2");
        ctx.topologyRegistry().addType("3");
        ctx.topologyRegistry().addType("4");
        ctx.topologyRegistry().addType("5");
        ctx.topologyRegistry().addType("6");
        ctx.topologyRegistry().addType("T1");
        ctx.topologyRegistry().addType("T2");
        ctx.topologyRegistry().addType("T3");
        ctx.topologyRegistry().addType("T4");
        ctx.particleTypes().add("p1", 1.);
        ctx.particleTypes().add("p2", 1.);
        ctx.particleTypes().add("p3", 1.);
        ctx.particleTypes().add("p4", 1.);

        readdy::model::top::reactions::STRParser parser (ctx.topologyRegistry());

        {
            auto r = parser.parse("topology-topology fusion type1: T1 (p1) + T2 (p2) -> T3 (p3--p4) [self=true]", 10., 11.);
            REQUIRE(r.mode() == readdy::model::top::reactions::STRMode::TT_FUSION_ALLOW_SELF);
            REQUIRE(r.top_type1() == ctx.topologyRegistry().idOf("T1"));
            REQUIRE(r.top_type2() == ctx.topologyRegistry().idOf("T2"));
            REQUIRE(r.type1() == ctx.particleTypes().idOf("p1"));
            REQUIRE(r.type2() == ctx.particleTypes().idOf("p2"));
            REQUIRE(r.top_type_to1() == ctx.topologyRegistry().idOf("T3"));
            REQUIRE(r.top_type_to2() == readdy::EmptyTopologyId);
            REQUIRE(r.type_to1() == ctx.particleTypes().idOf("p3"));
            REQUIRE(r.type_to2() == ctx.particleTypes().idOf("p4"));
            REQUIRE(r.name() == "topology-topology fusion type1");
            REQUIRE(r.rate() == 10.);
            REQUIRE(r.radius() == 11.);
	    REQUIRE(r.min_graph_distance() == 0);	    
        }
	{
	  auto r = parser.parse("topology-topology fusion type1: T1 (p1) + T2 (p2) -> T3 (p3--p4) [self=true, distance>6]", 10., 11.);
	  REQUIRE(r.mode() == readdy::model::top::reactions::STRMode::TT_FUSION_ALLOW_SELF);
	  REQUIRE(r.top_type1() == ctx.topologyRegistry().idOf("T1"));
	  REQUIRE(r.top_type2() == ctx.topologyRegistry().idOf("T2"));
	  REQUIRE(r.type1() == ctx.particleTypes().idOf("p1"));
	  REQUIRE(r.type2() == ctx.particleTypes().idOf("p2"));
	  REQUIRE(r.top_type_to1() == ctx.topologyRegistry().idOf("T3"));
	  REQUIRE(r.top_type_to2() == readdy::EmptyTopologyId);
	  REQUIRE(r.type_to1() == ctx.particleTypes().idOf("p3"));
	  REQUIRE(r.type_to2() == ctx.particleTypes().idOf("p4"));
	  REQUIRE(r.name() == "topology-topology fusion type1");
	  REQUIRE(r.rate() == 10.);
	  REQUIRE(r.radius() == 11.);
	  REQUIRE(r.min_graph_distance() == 6);
	}
	{
	  // distance option without self=true option should throw an error
	  try {
	    auto r = parser.parse("topology-topology fusion type1: T1 (p1) + T2 (p2) -> T3 (p3--p4) [distance>6]", 10., 11.);
	    REQUIRE( false );
	  } catch (std::runtime_error const &e) {
	    REQUIRE( true );
	  }
	}
        {
            auto r = parser.parse("topology-topology fusion type2: T1 (p1) + T2 (p2) -> T3 (p3--p4)", 10., 11.);
            REQUIRE(r.mode() == readdy::model::top::reactions::STRMode::TT_FUSION);
            REQUIRE(r.top_type1() == ctx.topologyRegistry().idOf("T1"));
            REQUIRE(r.top_type2() == ctx.topologyRegistry().idOf("T2"));
            REQUIRE(r.type1() == ctx.particleTypes().idOf("p1"));
            REQUIRE(r.type2() == ctx.particleTypes().idOf("p2"));
            REQUIRE(r.top_type_to1() == ctx.topologyRegistry().idOf("T3"));
            REQUIRE(r.top_type_to2() == readdy::EmptyTopologyId);
            REQUIRE(r.type_to1() == ctx.particleTypes().idOf("p3"));
            REQUIRE(r.type_to2() == ctx.particleTypes().idOf("p4"));
            REQUIRE(r.name() == "topology-topology fusion type2");
            REQUIRE(r.rate() == 10.);
            REQUIRE(r.radius() == 11.);
        }
        {
            auto r = parser.parse("topology-topology enzymatic type: T1 (p1) + T2 (p2) -> T3 (p3) + T4 (p4)", 10., 11.);
            REQUIRE(r.mode() == readdy::model::top::reactions::STRMode::TT_ENZYMATIC);
            REQUIRE(r.top_type1() == ctx.topologyRegistry().idOf("T1"));
            REQUIRE(r.top_type2() == ctx.topologyRegistry().idOf("T2"));
            REQUIRE(r.type1() == ctx.particleTypes().idOf("p1"));
            REQUIRE(r.type2() == ctx.particleTypes().idOf("p2"));
            REQUIRE(r.top_type_to1() == ctx.topologyRegistry().idOf("T3"));
            REQUIRE(r.top_type_to2() == ctx.topologyRegistry().idOf("T4"));
            REQUIRE(r.type_to1() == ctx.particleTypes().idOf("p3"));
            REQUIRE(r.type_to2() == ctx.particleTypes().idOf("p4"));
            REQUIRE(r.name() == "topology-topology enzymatic type");
            REQUIRE(r.rate() == 10.);
            REQUIRE(r.radius() == 11.);
        }
        {
            auto r = parser.parse("topology-particle fusion type: T1 (p1) + (p2) -> T2 (p3--p4)", 10., 11.);
            REQUIRE(r.mode() == readdy::model::top::reactions::STRMode::TP_FUSION);
            REQUIRE(r.top_type1() == ctx.topologyRegistry().idOf("T1"));
            REQUIRE(r.top_type2() == readdy::EmptyTopologyId);
            REQUIRE(r.type1() == ctx.particleTypes().idOf("p1"));
            REQUIRE(r.type2() == ctx.particleTypes().idOf("p2"));
            REQUIRE(r.top_type_to1() == ctx.topologyRegistry().idOf("T2"));
            REQUIRE(r.top_type_to2() == readdy::EmptyTopologyId);
            REQUIRE(r.type_to1() == ctx.particleTypes().idOf("p3"));
            REQUIRE(r.type_to2() == ctx.particleTypes().idOf("p4"));
            REQUIRE(r.name() == "topology-particle fusion type");
            REQUIRE(r.rate() == 10.);
            REQUIRE(r.radius() == 11.);
        }
        {
            auto r = parser.parse("topology-particle enzymatic type: T1 (p1) + (p2) -> T2 (p3) + (p4)", 10., 11.);
            REQUIRE(r.mode() == readdy::model::top::reactions::STRMode::TP_ENZYMATIC);
            REQUIRE(r.top_type1() == ctx.topologyRegistry().idOf("T1"));
            REQUIRE(r.top_type2() == readdy::EmptyTopologyId);
            REQUIRE(r.type1() == ctx.particleTypes().idOf("p1"));
            REQUIRE(r.type2() == ctx.particleTypes().idOf("p2"));
            REQUIRE(r.top_type_to1() == ctx.topologyRegistry().idOf("T2"));
            REQUIRE(r.top_type_to2() == readdy::EmptyTopologyId);
            REQUIRE(r.type_to1() == ctx.particleTypes().idOf("p3"));
            REQUIRE(r.type_to2() == ctx.particleTypes().idOf("p4"));
            REQUIRE(r.name() == "topology-particle enzymatic type");
            REQUIRE(r.rate() == 10.);
            REQUIRE(r.radius() == 11.);
        }
    }

    SECTION("Attach topologies") {
        ctx.periodicBoundaryConditions() = {{true, true, true}};
        ctx.topologyRegistry().addType("TA");
        ctx.boxSize() = {{15, 15, 15}};
        ctx.particleTypes().add("middle", 0., model::particleflavor::TOPOLOGY);
        ctx.particleTypes().add("end", 0., model::particleflavor::TOPOLOGY);
        ctx.topologyRegistry().configureBondPotential("middle", "middle", {.00000001, 1});
        ctx.topologyRegistry().configureBondPotential("middle", "end", {.00000001, 1});
        ctx.topologyRegistry().configureBondPotential("end", "end", {.00000001, 1});
        // register attach reaction that transforms (end, A) -> (middle, end)
        ctx.topologyRegistry().addSpatialReaction("merge: TA (end) + TA (end) -> TA (middle--middle)", 1e3, 1.5);

        Simulation simulation {create<TestType>(), ctx};

        for(int i = 0; i < 3; ++i) {
            auto top = simulation.addTopology("TA", {simulation.createTopologyParticle("end", {0. - 4. + 3*i, 0., 0.}),
                                                     simulation.createTopologyParticle("middle", {0. - 3. + 3*i, 0., 0.}),
                                                     simulation.createTopologyParticle("end", {0. - 2. + 3*i, 0., 0.})});
            {
                auto it = top->graph().vertices().begin();
                auto it2 = std::next(top->graph().vertices().begin());
                top->graph().addEdge(it, it2);
                ++it;
                ++it2;
                top->graph().addEdge(it, it2);
            }
        }

        

        REQUIRE(simulation.currentTopologies().size() == 3);
        simulation.run(6, 1.);

        const auto& type_registry = simulation.context().particleTypes();

        REQUIRE(simulation.context().topologyRegistry().isSpatialReactionType("end"));
        REQUIRE_FALSE(simulation.context().topologyRegistry().isSpatialReactionType("middle"));
        REQUIRE(simulation.context().calculateMaxCutoff() == 1.5);

        REQUIRE(simulation.currentTopologies().size() == 1);
        auto chainTop = simulation.currentTopologies().at(0);
        REQUIRE(chainTop->getNParticles() == 3 /*original topology particles*/ + 6 /*attached particles*/);

        auto top_particles = simulation.stateModel().getParticlesForTopology(*chainTop);

        bool foundEndVertex {false};
        // check that graph is indeed linear
        for(std::size_t idx = 0; idx < chainTop->graph().vertices().size() && !foundEndVertex; ++idx) {
            auto prev_neighbor = std::next(chainTop->graph().vertices().begin(), idx);
            const auto& v_end = *prev_neighbor;
            if(v_end.particleType() == type_registry.idOf("end")) {
                foundEndVertex = true;

                REQUIRE(v_end.neighbors().size() == 1);

                REQUIRE(top_particles.at(idx).type() == type_registry.idOf("end"));

                using flouble = fp::FloatingPoint<scalar>;
                flouble x_end (top_particles.at(idx).pos().x);
                flouble y_end (top_particles.at(idx).pos().y);
                flouble z_end (top_particles.at(idx).pos().z);
                // the end particle of our topology sausage should be either at x=4 or x=-4
                REQUIRE((x_end.AlmostEquals(flouble{4.}) || x_end.AlmostEquals(flouble{-4.})));
                REQUIRE(y_end.AlmostEquals(flouble{0.})); // no diffusion going on
                REQUIRE(z_end.AlmostEquals(flouble{0.})); // no diffusion going on

                auto factor = x_end.AlmostEquals(flouble{4.}) ? 1. : -1.;

                // walk along topology sausage, check end particles are always at +-4, the other ones are of type middle
                auto next_neighbor = v_end.neighbors().at(0);
                std::size_t i = 0;
                while(next_neighbor->particleType() != type_registry.idOf("end") && i < 20 /*no endless loop*/) {
                    auto next_idx = std::distance(chainTop->graph().vertices().begin(), next_neighbor);
                    const auto& next_particle = top_particles.at(static_cast<std::size_t>(next_idx));
                    auto predicted_pos = factor*4. - factor*(i+1)*1.;
                    auto actual_pos = next_particle.pos().x;
                    REQUIRE((flouble(actual_pos).AlmostEquals(flouble(predicted_pos))));
                    REQUIRE((flouble(next_particle.pos().y)).AlmostEquals(flouble(0.)));
                    REQUIRE((flouble(next_particle.pos().z)).AlmostEquals(flouble(0.)));
                    if(next_neighbor->particleType() == type_registry.idOf("middle")) {
                        REQUIRE(next_neighbor->neighbors().size() == 2);
                        if(next_neighbor->neighbors().at(0) == prev_neighbor) {
                            prev_neighbor = next_neighbor;
                            next_neighbor = next_neighbor->neighbors().at(1);
                        } else {
                            prev_neighbor = next_neighbor;
                            next_neighbor = next_neighbor->neighbors().at(0);
                        }

                    } else {
                        REQUIRE(next_neighbor->neighbors().size() == 1);
                    }

                    ++i;
                }

                REQUIRE(i == 3+6-1-1);
            }
        }
        REQUIRE(foundEndVertex);
    }
}
