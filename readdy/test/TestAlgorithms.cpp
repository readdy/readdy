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
 * @file TestAlgorithms.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 3/1/18
 */


#include <unordered_set>
#include <readdy/api/SimulationLoop.h>
#include <readdy/api/Simulation.h>
#include "gtest/gtest.h"
#include "readdy/common/algorithm.h"


namespace {

using namespace readdy;

struct Event {
    int i;
    scalar rate;
    scalar cumulativeRate;
};

TEST(TestLoop, Loop) {
Simulation sim ("SingleCPU");
sim.context().topologyRegistry().addType("Polymer");
sim.context().particleTypes().addTopologyType("A", 1);
sim.context().particleTypes().add("B", 1.);
sim.context().topologyRegistry().addSpatialReaction("Attach: Polymer(A) + (B) -> Polymer(A--A)", 1., 5.);
sim.context().reactions().add("myfus: B +(5) B -> B", 10);
sim.context().reactions().add("myconv: B -> B", 10);

sim.addParticle("B", 0., 0., 0.);
sim.addTopology("Polymer", {sim.createTopologyParticle("A", {0., 0., 0.})});
sim.context().topologyRegistry().configureBondPotential("A", "A", {10., 10.});

auto loop = sim.createLoop(1.5);

loop.addCallback([&](const auto t) {
    // triggers evaluation of rate functions for each topology
    for (auto topology : sim.stateModel().getTopologies()) {
        topology->updateReactionRates(
                sim.context().topologyRegistry().structuralReactionsOf(topology->type())
        );
    }
    // change spatial topology reaction rate
    sim.context().topologyRegistry().spatialReactionByName("Attach").rate() = 5000;
    // change unimolecular reaction rate
    sim.context().reactions().order1ByName("myconv")->rate() = 100;
    // change bimolecular reaction rate
    sim.context().reactions().order2ByName("myfus")->rate() = 100;
});

loop.run(500);
}

TEST(TestAlgorithms, EvaluateAll) {
    auto n = 1000U;
    std::vector<Event> events (n);
    for(auto i = 0U; i < n; ++i) {
        events.at(i).i = i;
        events.at(i).rate = i;
    }

    auto shouldEval = [](const Event &event) { return true; };
    auto depending = [](const Event &e1, const Event &e2) { return e1.rate == e2.rate; };

    std::unordered_set<int> set;
    auto eval = [&set](const Event &event) {
        set.insert(event.i);
    };

    algo::performEvents(events, shouldEval, depending, eval);

    for(auto i = 0; i < n; ++i) {
        ASSERT_TRUE(set.find(i) != set.end());
    }
}

TEST(TestAlgorithms, EvaluateEven) {
    auto n = 1000U;
    std::vector<Event> events (n);
    for(auto i = 0U; i < n; ++i) {
        events.at(i).i = i;
        events.at(i).rate = i/2;
    }

    auto shouldEval = [](const Event &event) { return true; };
    auto depending = [](const Event &e1, const Event &e2) { return e1.rate == e2.rate; };

    std::unordered_set<scalar> set;
    auto eval = [&set, &events](const Event &event) {
        {
            // check that the events are all moved to the end that were already taken care of
            for(const auto rate : set) {
                auto it1 = std::find_if(events.end() - 2*set.size(), events.end(), [rate](const Event &event) {
                    return event.rate == rate;
                });
                ASSERT_NE(it1, events.end());
                auto it2 = std::find_if(it1+1, events.end(), [rate](const Event &event) {
                    return event.rate == rate;
                });
                ASSERT_NE(it2, events.end());
            }
        }

        set.insert(event.rate);
    };

    algo::performEvents(events, shouldEval, depending, eval);

    for(auto i = 0; i < n/2; ++i) {
        ASSERT_TRUE(set.find(i) != set.end());
    }
}

}
