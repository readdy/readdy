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

class TestTopologyReactionsExternal : public KernelTest {
protected:
    void SetUp() override {
        auto &ctx = kernel->getKernelContext();
        ctx.particle_types().add("Topology A", 1.0, 1.0, readdy::model::particleflavor::TOPOLOGY);
        ctx.particle_types().add("Topology B", 1.0, 1.0, readdy::model::particleflavor::TOPOLOGY);
        ctx.particle_types().add("Topology Invalid Type", 1.0, 1.0, readdy::model::particleflavor::TOPOLOGY);
        ctx.particle_types().add("A", 1.0, 1.0, readdy::model::particleflavor::NORMAL);
        ctx.particle_types().add("B", 1.0, 1.0, readdy::model::particleflavor::NORMAL);

        ctx.configureTopologyBondPotential("Topology A", "Topology A", {10, 10});
        ctx.configureTopologyBondPotential("Topology A", "Topology B", {10, 10});
        ctx.configureTopologyBondPotential("Topology B", "Topology B", {10, 10});

        ctx.setBoxSize(10, 10, 10);
    }
};

namespace {

TEST_P(TestTopologyReactionsExternal, TestTopologyEnzymaticReaction) {
    using namespace readdy;
    auto &ctx = kernel->getKernelContext();
    model::TopologyParticle x_0{c_::zero, c_::zero, c_::zero, ctx.particle_types().id_of("Topology A")};
    kernel->getKernelStateModel().addTopology({x_0});
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
    auto toplogy = kernel->getKernelStateModel().addTopology({x_0});
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
    sim.setKernel(kernel->getName());
    sim.setPeriodicBoundary({{true, true, true}});
    sim.setBoxSize(100, 100, 100);
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
    auto toplogy = sim.addTopology(topologyParticles);
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

    model::top::reactions::TopologyReaction r {[aId](model::top::GraphTopology& top) {
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

    toplogy->addReaction(std::move(r));

    sim.runScheme<api::ReaDDyScheme>(true).evaluateTopologyReactions().configureAndRun(35, 1.);

    log::debug("got n topologies: {}", sim.currentTopologies().size());
    const auto kernel = sim.getSelectedKernel();
    for(auto top : sim.currentTopologies()) {
        for(const auto p : top->getParticles()) {
            ASSERT_EQ(kernel->getKernelStateModel().getTopologyForParticle(p), top);
        }
    }

}

INSTANTIATE_TEST_CASE_P(Kernels, TestTopologyReactionsExternal,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}
