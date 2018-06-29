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
#include <readdy/model/reactions/Reaction.h>
#include <readdy/testing/NOOPPotential.h>

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
    context.particleTypes().add("A", 1.);
    auto &compartments = context.compartments();
    auto id = compartments.addPlane({{"A", "A"}}, "Plane", {1., 0., 0.}, 10., true);
    ASSERT_EQ(compartments.get().size(), 1);
    ASSERT_EQ(compartments.get().begin()->get()->getId(), id);
}

TEST(Context, Topologies) {
    Context context;
    context.particleTypes().addTopologyType("A", 1.);
    context.particleTypes().addTopologyType("B", 1.);
    context.particleTypes().addTopologyType("C", 1.);
    context.particleTypes().add("P", 1.);
    context.particleTypes().add("Q", 1.);
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
        ASSERT_EQ(reaction.type1(), context.particleTypes().idOf("A"));
        ASSERT_EQ(reaction.type2(), context.particleTypes().idOf("P"));
        ASSERT_EQ(reaction.top_type1(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type_to1(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.type_to1(), context.particleTypes().idOf("B"));
        ASSERT_EQ(reaction.type_to2(), context.particleTypes().idOf("Q"));
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
        ASSERT_EQ(reaction.type1(), context.particleTypes().idOf("A"));
        ASSERT_EQ(reaction.type2(), context.particleTypes().idOf("B"));
        ASSERT_EQ(reaction.top_type1(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type2(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type_to1(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.top_type_to2(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.type_to1(), context.particleTypes().idOf("B"));
        ASSERT_EQ(reaction.type_to2(), context.particleTypes().idOf("A"));
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
        ASSERT_EQ(reaction.type1(), context.particleTypes().idOf("A"));
        ASSERT_EQ(reaction.type2(), context.particleTypes().idOf("P"));
        ASSERT_EQ(reaction.top_type1(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type_to1(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.type_to1(), context.particleTypes().idOf("B"));
        ASSERT_EQ(reaction.type_to2(), context.particleTypes().idOf("B"));
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
        ASSERT_EQ(reaction.type1(), context.particleTypes().idOf("A"));
        ASSERT_EQ(reaction.type2(), context.particleTypes().idOf("B"));
        ASSERT_EQ(reaction.top_type1(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type2(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type_to1(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.type_to1(), context.particleTypes().idOf("B"));
        ASSERT_EQ(reaction.type_to2(), context.particleTypes().idOf("A"));
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
        ASSERT_EQ(reaction.type1(), context.particleTypes().idOf("A"));
        ASSERT_EQ(reaction.type2(), context.particleTypes().idOf("B"));
        ASSERT_EQ(reaction.top_type1(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type2(), topologies.idOf("T"));
        ASSERT_EQ(reaction.top_type_to1(), topologies.idOf("TT"));
        ASSERT_EQ(reaction.type_to1(), context.particleTypes().idOf("B"));
        ASSERT_EQ(reaction.type_to2(), context.particleTypes().idOf("A"));
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


TEST(Context, ReactionDescriptorAddReactions) {
    const auto prepareCtx = [](Context &ctx, const std::string& descriptor, readdy::scalar rate){
        ctx.particleTypes().add("A", 1.);
        ctx.particleTypes().add("B", 1.);
        ctx.particleTypes().add("C", 1.);
        ctx.reactions().add(descriptor, rate);
    };
    {
        Context ctx;
        auto decay = "mydecay:A->";
        prepareCtx(ctx, decay, 2.);
        const auto &r = ctx.reactions().order1ByName("mydecay");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->type(), readdy::model::reactions::ReactionType::Decay);
        EXPECT_EQ(r->nEducts(), 1);
        EXPECT_EQ(r->nProducts(), 0);
        EXPECT_EQ(r->educts()[0], ctx.particleTypes().idOf("A"));
        EXPECT_EQ(r->rate(), 2.);
    }
    {
        Context ctx;
        auto conversion = "myconv: A -> B";
        prepareCtx(ctx, conversion, 3.);
        const auto &r = ctx.reactions().order1ByName("myconv");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->type(), readdy::model::reactions::ReactionType::Conversion);
        EXPECT_EQ(r->nEducts(), 1);
        EXPECT_EQ(r->nProducts(), 1);
        EXPECT_EQ(r->educts()[0], ctx.particleTypes().idOf("A"));
        EXPECT_EQ(r->products()[0], ctx.particleTypes().idOf("B"));
        EXPECT_EQ(r->rate(), 3.);
    }
    {
        Context ctx;
        auto fusion = "myfus: B +(1.2) B -> C [0.5, 0.5]";
        prepareCtx(ctx, fusion, 4.);
        const auto &r = ctx.reactions().order2ByName("myfus");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->type(), readdy::model::reactions::ReactionType::Fusion);
        EXPECT_EQ(r->nEducts(), 2);
        EXPECT_EQ(r->nProducts(), 1);
        EXPECT_EQ(r->educts()[0], ctx.particleTypes().idOf("B"));
        EXPECT_EQ(r->educts()[1], ctx.particleTypes().idOf("B"));
        EXPECT_EQ(r->products()[0], ctx.particleTypes().idOf("C"));
        EXPECT_EQ(r->eductDistance(), 1.2);
        EXPECT_EQ(r->weight1(), 0.5);
        EXPECT_EQ(r->weight2(), 0.5);
        EXPECT_EQ(r->rate(), 4.);
    }
    {
        Context ctx;
        auto fission = "myfiss: B -> C +(3.0) B [0.1, 0.9]";
        prepareCtx(ctx, fission, 5.);
        const auto &r = ctx.reactions().order1ByName("myfiss");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->type(), readdy::model::reactions::ReactionType::Fission);
        EXPECT_EQ(r->nEducts(), 1);
        EXPECT_EQ(r->nProducts(), 2);
        EXPECT_EQ(r->educts()[0], ctx.particleTypes().idOf("B"));
        EXPECT_EQ(r->products()[0], ctx.particleTypes().idOf("C"));
        EXPECT_EQ(r->products()[1], ctx.particleTypes().idOf("B"));
        EXPECT_EQ(r->productDistance(), 3.0);
        EXPECT_EQ(r->weight1(), 0.1);
        EXPECT_EQ(r->weight2(), 0.9);
        EXPECT_EQ(r->rate(), 5.);
    }
    {
        Context ctx;
        auto enzymatic = "myenz:A +(1.5) C -> B + C";
        prepareCtx(ctx, enzymatic, 6.);
        const auto &r = ctx.reactions().order2ByName("myenz");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->type(), readdy::model::reactions::ReactionType::Enzymatic);
        EXPECT_EQ(r->nEducts(), 2);
        EXPECT_EQ(r->nProducts(), 2);
        EXPECT_EQ(r->educts()[0], ctx.particleTypes().idOf("A"));
        EXPECT_EQ(r->educts()[1], ctx.particleTypes().idOf("C"));
        EXPECT_EQ(r->products()[0], ctx.particleTypes().idOf("B"));
        EXPECT_EQ(r->products()[1], ctx.particleTypes().idOf("C"));
        EXPECT_EQ(r->eductDistance(), 1.5);
        EXPECT_EQ(r->rate(), 6.);
    }
    {
        Context ctx;
        auto enzymatic = "eat:B +(1) A -> A + A";
        prepareCtx(ctx, enzymatic, 7.);
        const auto &r = ctx.reactions().order2ByName("eat");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->type(), readdy::model::reactions::ReactionType::Enzymatic);
        EXPECT_EQ(r->nEducts(), 2);
        EXPECT_EQ(r->nProducts(), 2);
        EXPECT_EQ(r->educts()[0], ctx.particleTypes().idOf("B"));
        EXPECT_EQ(r->educts()[1], ctx.particleTypes().idOf("A"));
        EXPECT_EQ(r->products()[0], ctx.particleTypes().idOf("A"));
        EXPECT_EQ(r->products()[1], ctx.particleTypes().idOf("A"));
        EXPECT_EQ(r->eductDistance(), 1.);
        EXPECT_EQ(r->rate(), 7.);
    }
}

TEST(Context, ReactionDescriptorInvalidInputs) {
    Context ctx;
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    std::vector<std::string> inv = {"myinvalid: + A -> B", "noarrow: A B", " : Noname ->", "weights: A + A -> A [0.1, ]", "blub: A (3)+ A -> B", "eat: A +(1) B -> A + A"};
    for (const auto &i : inv) {
        EXPECT_ANY_THROW(ctx.reactions().add(i, 42.));
    }
    EXPECT_EQ(ctx.reactions().nOrder1(), 0);
    EXPECT_EQ(ctx.reactions().nOrder2(), 0);
}

TEST(Context, ReactionNameExists) {
    Context ctx;
    ctx.particleTypes().add("A", 1.);
    ctx.reactions().add("bla: A->", 1.);
    EXPECT_ANY_THROW(ctx.reactions().add("bla: A->", 1.));
}

TEST(Context, ReactionNameAndId) {
    Context ctx;
    ctx.particleTypes().add("A", 1.);
    ctx.reactions().add("foo: A->", 1.);
    ctx.reactions().add("bla: A+(1)A->A", 1.);
    const auto idFoo = ctx.reactions().idOf("foo");
    const auto idBla = ctx.reactions().idOf("bla");
    EXPECT_GE(idFoo, 0);
    EXPECT_GE(idBla, 0);
    EXPECT_EQ(ctx.reactions().nameOf(idFoo), "foo");
    EXPECT_EQ(ctx.reactions().nameOf(idBla), "bla");
}

TEST(Context, GetAllReactions) {
    Context ctx;
    ctx.particleTypes().add("A", 1.);
    ctx.reactions().add("foo: A->", 1.);
    ctx.reactions().add("bla: A+(1)A->A", 1.);
    ctx.reactions().addConversion("conv1", "A", "A", 1.);
    ctx.reactions().addConversion("conv2", "A", "A", 1.);
    ctx.reactions().addFusion("fusion", "A","A", "A", 1., 1.);
    const auto &o1flat = ctx.reactions().order1Flat();
    const auto &o2flat = ctx.reactions().order2Flat();
    EXPECT_EQ(o1flat.size() + o2flat.size(), 5);
}

TEST(Context, GetReactionById) {
    Context ctx;
    ctx.particleTypes().add("A", 1.);
    ctx.reactions().add("foo: A->", 1.);
    ctx.reactions().add("bla: A+(1)A->A", 1.);
    auto idFoo = ctx.reactions().idOf("foo");
    auto idBla = ctx.reactions().idOf("bla");
    auto foo = ctx.reactions().byId(idFoo);
    auto bla = ctx.reactions().byId(idBla);
    EXPECT_EQ(foo->name(), "foo");
    EXPECT_EQ(bla->name(), "bla");
}

TEST(Context, PotentialOrder2Map) {
    Context ctx;
    ctx.particleTypes().add("a", 1.);
    ctx.particleTypes().add("b", 1.);
    auto noop = std::make_unique<readdy::testing::NOOPPotentialOrder2>(
            ctx.particleTypes()("a"), ctx.particleTypes()("b"));
    ctx.potentials().addUserDefined(noop.get());
    auto noop2 = std::make_unique<readdy::testing::NOOPPotentialOrder2>(
            ctx.particleTypes()("b"), ctx.particleTypes()("a"));
    ctx.potentials().addUserDefined(noop2.get());

    ctx.potentials().addScreenedElectrostatics("a", "b", 1., 1., 1., 1., 3, .1);
    {
        const auto &vector = ctx.potentials().potentialsOf("b", "a");
        EXPECT_EQ(vector.size(), 3);
        EXPECT_EQ(vector[0], noop.get());
        EXPECT_EQ(vector[1], noop2.get());
    }
    {
        const auto &vector = ctx.potentials().potentialsOf("a", "b");
        EXPECT_EQ(vector.size(), 3);
        EXPECT_EQ(vector[0], noop.get());
        EXPECT_EQ(vector[1], noop2.get());
    }
}

TEST(TestKernelContext, Copyable) {
    Context context;
    Context copy(context);
}


}
