/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * @file TestContext.cpp
 * @brief Tests for the context class.
 * @author clonker
 * @date 6/25/18
 */

#include "gtest/gtest.h"
#include <readdy/model/Context.h>
#include <gmock/gmock-matchers.h>

namespace {

using Context = readdy::model::Context;

TEST(Context, KBT) {
    Context context;
    context.kBT() = 42;
    ASSERT_EQ(context.kBT(), 42);
}

TEST(Context, Box) {
    Context context;
    context.boxSize() = {{10, 20, 30}};
    context.periodicBoundaryConditions() = {{true, false, true}};
    ASSERT_EQ(context.boxVolume(), 10*20*30);
    {
        std::array<readdy::scalar, 3> box{{10, 20, 30}};
        ASSERT_EQ(context.boxSize(), box);
    }
    {
        std::array<bool, 3> pbc {{true, false, true}};
        ASSERT_EQ(context.periodicBoundaryConditions(), pbc);
    }
    {
        auto bbox = context.getBoxBoundingVertices();
        ASSERT_EQ(std::get<0>(bbox), readdy::Vec3(-5, -10, -15));
        ASSERT_EQ(std::get<1>(bbox), readdy::Vec3(5, 10, 15));
    }
}

TEST(Context, RecordStuff) {
    Context context;
    context.recordReactionCounts() = true;
    context.recordReactionsWithPositions() = true;
    context.recordVirial() = true;
    ASSERT_TRUE(context.recordReactionCounts());
    ASSERT_TRUE(context.recordReactionsWithPositions());
    ASSERT_TRUE(context.recordVirial());
    context.recordReactionCounts() = false;
    context.recordReactionsWithPositions() = false;
    context.recordVirial() = false;
    ASSERT_FALSE(context.recordReactionCounts());
    ASSERT_FALSE(context.recordReactionsWithPositions());
    ASSERT_FALSE(context.recordVirial());
}

TEST(Context, Compartments) {
    Context context;
    context.particle_types().add("A", 1.);
    auto &compartments = context.compartments();
    auto id = compartments.addPlane({{"A", "A"}}, "Plane", {1., 0., 0.}, 10., true);
    ASSERT_EQ(compartments.get().size(), 1);
    ASSERT_EQ(compartments.get().begin()->get()->getId(), id);
}

TEST(Context, Topologies) {
    Context context;
    context.particle_types().addTopologyType("A", 1.);
    context.particle_types().addTopologyType("B", 1.);
    context.particle_types().addTopologyType("C", 1.);
    context.particle_types().add("P", 1.);
    context.particle_types().add("Q", 1.);
    auto &topologies = context.topologyRegistry();

    topologies.addType("T");
    topologies.addType("TT");
    ASSERT_EQ(topologies.types().size(), 2);
    ASSERT_EQ(topologies.types()[0].name, "T");
    ASSERT_EQ(topologies.types()[1].name, "TT");

    topologies.addSpatialReaction("MyFusionReaction: T(A)+T(B) -> TT(B--A)", 10, 20);
    topologies.addSpatialReaction("MySelfFusionReaction: T(A)+T(B) -> TT(B--A) [self=true]", 10, 20);
    topologies.addSpatialReaction("MyTPFusion: T(A)+(P) -> TT(B--B)", 10, 20);
    topologies.addSpatialReaction("MyEnzymaticTT: T(A)+T(B) -> TT(B) + TT(A)", 10, 20);
    topologies.addSpatialReaction("MyEnzymaticTP: T(A) + (P) -> TT(B) + (Q)", 10, 20);

    ASSERT_EQ(topologies.spatialReactionRegistry().size(), 2);
    ASSERT_TRUE(topologies.isSpatialReactionType("A"));
    ASSERT_FALSE(topologies.isSpatialReactionType("C"));

    {
        // MyEnzymaticTP
        auto reaction = topologies.spatialReactionByName("MyEnzymaticTP");
        ASSERT_EQ(reaction.name(), "MyEnzymaticTP");
        ASSERT_EQ(reaction.type1(), context.particle_types().idOf("A"));
        ASSERT_EQ(reaction.type2(), context.particle_types().idOf("P"));
        ASSERT_EQ(reaction.top_type1(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type_to1(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.type_to1(), context.particle_types().idOf("B"));
        ASSERT_EQ(reaction.type_to2(), context.particle_types().idOf("Q"));
        ASSERT_TRUE(reaction.is_topology_particle_reaction());
        ASSERT_FALSE(reaction.is_topology_topology_reaction());
        ASSERT_TRUE(reaction.is_enzymatic());
        ASSERT_FALSE(reaction.is_fusion());
        ASSERT_EQ(reaction.rate(), 10);
        ASSERT_EQ(reaction.radius(), 20);
        ASSERT_FALSE(reaction.allow_self_connection());
        ASSERT_EQ(reaction.mode(), readdy::model::top::reactions::STRMode::TP_ENZYMATIC);
    }

    {
        // MyEnzymaticTT
        auto reaction = topologies.spatialReactionByName("MyEnzymaticTT");
        ASSERT_EQ(reaction.name(), "MyEnzymaticTT");
        ASSERT_EQ(reaction.type1(), context.particle_types().idOf("A"));
        ASSERT_EQ(reaction.type2(), context.particle_types().idOf("B"));
        ASSERT_EQ(reaction.top_type1(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type2(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type_to1(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.top_type_to2(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.type_to1(), context.particle_types().idOf("B"));
        ASSERT_EQ(reaction.type_to2(), context.particle_types().idOf("A"));
        ASSERT_FALSE(reaction.is_topology_particle_reaction());
        ASSERT_TRUE(reaction.is_topology_topology_reaction());
        ASSERT_TRUE(reaction.is_enzymatic());
        ASSERT_FALSE(reaction.is_fusion());
        ASSERT_EQ(reaction.rate(), 10);
        ASSERT_EQ(reaction.radius(), 20);
        ASSERT_FALSE(reaction.allow_self_connection());
        ASSERT_EQ(reaction.mode(), readdy::model::top::reactions::STRMode::TT_ENZYMATIC);
    }

    {
        // MyTPFusion
        auto reaction = topologies.spatialReactionByName("MyTPFusion");
        ASSERT_EQ(reaction.name(), "MyTPFusion");
        ASSERT_EQ(reaction.type1(), context.particle_types().idOf("A"));
        ASSERT_EQ(reaction.type2(), context.particle_types().idOf("P"));
        ASSERT_EQ(reaction.top_type1(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type_to1(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.type_to1(), context.particle_types().idOf("B"));
        ASSERT_EQ(reaction.type_to2(), context.particle_types().idOf("B"));
        ASSERT_TRUE(reaction.is_topology_particle_reaction());
        ASSERT_FALSE(reaction.is_topology_topology_reaction());
        ASSERT_FALSE(reaction.is_enzymatic());
        ASSERT_TRUE(reaction.is_fusion());
        ASSERT_EQ(reaction.rate(), 10);
        ASSERT_EQ(reaction.radius(), 20);
        ASSERT_FALSE(reaction.allow_self_connection());
        ASSERT_EQ(reaction.mode(), readdy::model::top::reactions::STRMode::TP_FUSION);
    }

    {
        // MyFusionReaction
        auto reaction = topologies.spatialReactionByName("MyFusionReaction");
        ASSERT_EQ(reaction.name(), "MyFusionReaction");
        ASSERT_EQ(reaction.type1(), context.particle_types().idOf("A"));
        ASSERT_EQ(reaction.type2(), context.particle_types().idOf("B"));
        ASSERT_EQ(reaction.top_type1(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type2(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type_to1(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.type_to1(), context.particle_types().idOf("B"));
        ASSERT_EQ(reaction.type_to2(), context.particle_types().idOf("A"));
        ASSERT_FALSE(reaction.is_topology_particle_reaction());
        ASSERT_TRUE(reaction.is_topology_topology_reaction());
        ASSERT_FALSE(reaction.is_enzymatic());
        ASSERT_TRUE(reaction.is_fusion());
        ASSERT_EQ(reaction.rate(), 10);
        ASSERT_EQ(reaction.radius(), 20);
        ASSERT_FALSE(reaction.allow_self_connection());
        ASSERT_EQ(reaction.mode(), readdy::model::top::reactions::STRMode::TT_FUSION);
    }
    {
        // MySelfFusionReaction
        auto reaction = topologies.spatialReactionByName("MySelfFusionReaction");
        ASSERT_EQ(reaction.name(), "MySelfFusionReaction");
        ASSERT_EQ(reaction.type1(), context.particle_types().idOf("A"));
        ASSERT_EQ(reaction.type2(), context.particle_types().idOf("B"));
        ASSERT_EQ(reaction.top_type1(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type2(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type_to1(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.type_to1(), context.particle_types().idOf("B"));
        ASSERT_EQ(reaction.type_to2(), context.particle_types().idOf("A"));
        ASSERT_FALSE(reaction.is_topology_particle_reaction());
        ASSERT_TRUE(reaction.is_topology_topology_reaction());
        ASSERT_FALSE(reaction.is_enzymatic());
        ASSERT_TRUE(reaction.is_fusion());
        ASSERT_EQ(reaction.rate(), 10);
        ASSERT_EQ(reaction.radius(), 20);
        ASSERT_TRUE(reaction.allow_self_connection());
        ASSERT_EQ(reaction.mode(), readdy::model::top::reactions::STRMode::TT_FUSION_ALLOW_SELF);
    }
}

}
