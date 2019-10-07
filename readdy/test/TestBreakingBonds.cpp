/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
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
 * @file TestBreakingBonds.cpp
 * @brief Test implementation-independent execution of action that breaks bonds
 * @author chrisfroe
 * @date 27.05.19
 */

#include <catch2/catch.hpp>

#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

namespace m = readdy::model;
using namespace readdytesting::kernel;

TEMPLATE_TEST_CASE("Test breaking bonds.", "[breakbonds]", SingleCPU, CPU) {
    readdy::log::set_level(spdlog::level::warn);

    auto kernel = readdytesting::kernel::create<TestType>();
    auto &ctx = kernel->context();
    ctx.boxSize() = {10., 10., 10.};
    auto &types = ctx.particleTypes();
    auto &stateModel = kernel->stateModel();

    types.add("A", 1.0);
    types.add("B", 1.0);

    auto &topReg = ctx.topologyRegistry();
    topReg.addType("T");
    readdy::scalar timeStep = 1.;

    SECTION("Dimer") {
        readdy::api::Bond bond{1., 1., readdy::api::BondType::HARMONIC};
        topReg.configureBondPotential("A", "A", bond);

        std::vector<readdy::model::TopologyParticle> particles{
                {0., 0., 0., types.idOf("A")},
                {0., 0., 2., types.idOf("A")}
        };
        // distance of A and A is 2, equilibrium distance is 1, bond extension is 1, bond energy is 1.

        auto graphTop = stateModel.addTopology(topReg.idOf("T"), particles);
        auto &graph = graphTop->graph();
        graph.addEdgeBetweenParticles(0, 1);

        // We have one dimer
        auto topsBefore = stateModel.getTopologies();
        REQUIRE(topsBefore.size() == 1);
        REQUIRE(topsBefore.at(0)->type() == topReg.idOf("T"));
        REQUIRE(topsBefore.at(0)->getNParticles() == 2);

        SECTION("Low threshold, break") {
            readdy::model::actions::top::BreakConfig breakConfig;
            breakConfig.addBreakablePair(types.idOf("A"), types.idOf("A"), 0.9, 1e10);

            auto &&breakingBonds = kernel->actions().breakBonds(timeStep, breakConfig);
            breakingBonds->perform();

            // Bond broke, thus no topologies anymore
            auto topsAfter = stateModel.getTopologies();
            //readdy::log::warn("topology[0] has size {}", topsAfter.at(0)->getNParticles());
            REQUIRE(topsAfter.empty());
        }

        SECTION("High threshold, do not break") {
            readdy::model::actions::top::BreakConfig breakConfig;
            breakConfig.addBreakablePair(types.idOf("A"), types.idOf("A"), 1.1, 1e10);

            auto &&breakingBonds = kernel->actions().breakBonds(timeStep, breakConfig);
            breakingBonds->perform();

            // We still have one dimer
            auto topsAfter = stateModel.getTopologies();
            REQUIRE(topsAfter.size() == 1);
            REQUIRE(topsAfter.at(0)->type() == topReg.idOf("T"));
            REQUIRE(topsAfter.at(0)->getNParticles() == 2);
        }

        // particles are still at their position
        auto ps = stateModel.getParticles();
        REQUIRE(ps.size() == 2);
        REQUIRE(ps.at(0).type() == types.idOf("A"));
        REQUIRE(ps.at(1).type() == types.idOf("A"));
        if (readdy::testing::vec3eq(ps.at(0).pos(), readdy::Vec3(0., 0., 0.))) {
            REQUIRE(readdy::testing::vec3eq(ps.at(1).pos(), readdy::Vec3(0., 0., 2.)));
        } else {
            REQUIRE(readdy::testing::vec3eq(ps.at(0).pos(), readdy::Vec3(0., 0., 2.)));
            REQUIRE(readdy::testing::vec3eq(ps.at(1).pos(), readdy::Vec3(0., 0., 0.)));
        }

    }

    SECTION("Break due to large displacement with high rate, trimer") {

    }

    SECTION("Break due to large displacement with high rate, different types") {

    }

    SECTION("Break due to umbrella (box potential) pulling the bond apart (integration w/ forces and diffusion)") {

    }

}
