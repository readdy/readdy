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
 * @file TestContext.cpp
 * @brief Tests for the context class.
 * @author clonker
 * @date 6/25/18
 */

#include <string_view>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <readdy/model/Context.h>
#include <readdy/model/reactions/Reaction.h>
#include <readdy/testing/NOOPPotential.h>

using Context = readdy::model::Context;

TEST_CASE("Test context.", "[context]") {
    Context context;
    SECTION("kBT") {
        context.kBT() = 42;
        REQUIRE(context.kBT() == 42);
    }
    SECTION("Simulation box settings") {
        context.boxSize() = {{10, 20, 30}};
        context.periodicBoundaryConditions() = {{true, false, true}};
        REQUIRE(context.boxVolume() == 10*20*30);
        {
            std::array<readdy::scalar, 3> box{{10, 20, 30}};
            REQUIRE(context.boxSize() == box);
        }
        {
            std::array<bool, 3> pbc {{true, false, true}};
            REQUIRE(context.periodicBoundaryConditions() == pbc);
        }
        {
            auto bbox = context.getBoxBoundingVertices();
            REQUIRE(std::get<0>(bbox) == readdy::Vec3(-5, -10, -15));
            REQUIRE(std::get<1>(bbox) == readdy::Vec3(5, 10, 15));
        }
    }
    SECTION("Observable record switches") {
        context.recordReactionCounts() = true;
        context.recordReactionsWithPositions() = true;
        context.recordVirial() = true;
        REQUIRE(context.recordReactionCounts());
        REQUIRE(context.recordReactionsWithPositions());
        REQUIRE(context.recordVirial());
        context.recordReactionCounts() = false;
        context.recordReactionsWithPositions() = false;
        context.recordVirial() = false;
        REQUIRE_FALSE(context.recordReactionCounts());
        REQUIRE_FALSE(context.recordReactionsWithPositions());
        REQUIRE_FALSE(context.recordVirial());
    }
    SECTION("Compartments") {
        context.particleTypes().add("A", 1.);
        auto &compartments = context.compartments();
        auto id = compartments.addPlane({{"A", "A"}}, "Plane", {1., 0., 0.}, 10., true);
        REQUIRE(compartments.get().size() == 1);
        REQUIRE(compartments.get().begin()->get()->getId() == id);
    }

    SECTION("Topologies") {
        context.particleTypes().addTopologyType("A", 1.);
        context.particleTypes().addTopologyType("B", 1.);
        context.particleTypes().addTopologyType("C", 1.);
        context.particleTypes().add("P", 1.);
        context.particleTypes().add("Q", 1.);
        auto &topologies = context.topologyRegistry();

        topologies.addType("T");
        topologies.addType("TT");
        REQUIRE(topologies.types().size() == 2);
        REQUIRE(topologies.types()[0].name == "T");
        REQUIRE(topologies.types()[1].name == "TT");

        topologies.addSpatialReaction("MyFusionReaction: T(A)+T(B) -> TT(B--A)", 10, 20);
        topologies.addSpatialReaction("MySelfFusionReaction: T(A)+T(B) -> TT(B--A) [self=true]", 10, 20);
        topologies.addSpatialReaction("MyTPFusion: T(A)+(P) -> TT(B--B)", 10, 20);
        topologies.addSpatialReaction("MyEnzymaticTT: T(A)+T(B) -> TT(B) + TT(A)", 10, 20);
        topologies.addSpatialReaction("MyEnzymaticTP: T(A) + (P) -> TT(B) + (Q)", 10, 20);

        CHECK(topologies.spatialReactionRegistry().size() == 2);
        CHECK(topologies.isSpatialReactionType("A"));
        CHECK_FALSE(topologies.isSpatialReactionType("C"));

        SECTION("TP_Enzymatic") {
            auto reaction = topologies.spatialReactionByName("MyEnzymaticTP");
            CHECK(reaction.name() == "MyEnzymaticTP");
            CHECK(reaction.type1() == context.particleTypes().idOf("A"));
            CHECK(reaction.type2() == context.particleTypes().idOf("P"));
            CHECK(reaction.top_type1() == topologies.idOf("T"));
            CHECK(reaction.top_type_to1() == topologies.idOf("TT"));
            CHECK(reaction.type_to1() == context.particleTypes().idOf("B"));
            CHECK(reaction.type_to2() == context.particleTypes().idOf("Q"));
            CHECK(reaction.is_topology_particle_reaction());
            CHECK_FALSE(reaction.is_topology_topology_reaction());
            CHECK(reaction.is_enzymatic());
            CHECK_FALSE(reaction.is_fusion());
            CHECK(reaction.rate() == 10);
            CHECK(reaction.radius() == 20);
            CHECK_FALSE(reaction.allow_self_connection());
            CHECK(reaction.mode() == readdy::model::top::reactions::STRMode::TP_ENZYMATIC);
        }
        SECTION("TT_Enzymatic") {
            auto reaction = topologies.spatialReactionByName("MyEnzymaticTT");
            CHECK(reaction.name() == "MyEnzymaticTT");
            CHECK(reaction.type1() == context.particleTypes().idOf("A"));
            CHECK(reaction.type2() == context.particleTypes().idOf("B"));
            CHECK(reaction.top_type1() == topologies.idOf("T"));
            CHECK(reaction.top_type2() == topologies.idOf("T"));
            CHECK(reaction.top_type_to1() == topologies.idOf("TT"));
            CHECK(reaction.top_type_to2() == topologies.idOf("TT"));
            CHECK(reaction.type_to1() == context.particleTypes().idOf("B"));
            CHECK(reaction.type_to2() == context.particleTypes().idOf("A"));
            CHECK_FALSE(reaction.is_topology_particle_reaction());
            CHECK(reaction.is_topology_topology_reaction());
            CHECK(reaction.is_enzymatic());
            CHECK_FALSE(reaction.is_fusion());
            CHECK(reaction.rate() == 10);
            CHECK(reaction.radius() == 20);
            CHECK_FALSE(reaction.allow_self_connection());
            CHECK(reaction.mode() == readdy::model::top::reactions::STRMode::TT_ENZYMATIC);
        }
        SECTION("TP_Fusion") {
            auto reaction = topologies.spatialReactionByName("MyTPFusion");
            CHECK(reaction.name() == "MyTPFusion");
            CHECK(reaction.type1() == context.particleTypes().idOf("A"));
            CHECK(reaction.type2() == context.particleTypes().idOf("P"));
            CHECK(reaction.top_type1() == topologies.idOf("T"));
            CHECK(reaction.top_type_to1() == topologies.idOf("TT"));
            CHECK(reaction.type_to1() == context.particleTypes().idOf("B"));
            CHECK(reaction.type_to2() == context.particleTypes().idOf("B"));
            CHECK(reaction.is_topology_particle_reaction());
            CHECK_FALSE(reaction.is_topology_topology_reaction());
            CHECK_FALSE(reaction.is_enzymatic());
            CHECK(reaction.is_fusion());
            CHECK(reaction.rate() == 10);
            CHECK(reaction.radius() == 20);
            CHECK_FALSE(reaction.allow_self_connection());
            CHECK(reaction.mode() == readdy::model::top::reactions::STRMode::TP_FUSION);
        }
        SECTION("TT_Fusion") {
            auto reaction = topologies.spatialReactionByName("MyFusionReaction");
            CHECK(reaction.name() == "MyFusionReaction");
            CHECK(reaction.type1() == context.particleTypes().idOf("A"));
            CHECK(reaction.type2() == context.particleTypes().idOf("B"));
            CHECK(reaction.top_type1() == topologies.idOf("T"));
            CHECK(reaction.top_type2() == topologies.idOf("T"));
            CHECK(reaction.top_type_to1() == topologies.idOf("TT"));
            CHECK(reaction.type_to1() == context.particleTypes().idOf("B"));
            CHECK(reaction.type_to2() == context.particleTypes().idOf("A"));
            CHECK_FALSE(reaction.is_topology_particle_reaction());
            CHECK(reaction.is_topology_topology_reaction());
            CHECK_FALSE(reaction.is_enzymatic());
            CHECK(reaction.is_fusion());
            CHECK(reaction.rate() == 10);
            CHECK(reaction.radius() == 20);
            CHECK_FALSE(reaction.allow_self_connection());
            CHECK(reaction.mode() == readdy::model::top::reactions::STRMode::TT_FUSION);
        }
        SECTION("TT_Fusion_allow_self") {
            auto reaction = topologies.spatialReactionByName("MySelfFusionReaction");
            CHECK(reaction.name() == "MySelfFusionReaction");
            CHECK(reaction.type1() == context.particleTypes().idOf("A"));
            CHECK(reaction.type2() == context.particleTypes().idOf("B"));
            CHECK(reaction.top_type1() == topologies.idOf("T"));
            CHECK(reaction.top_type2() == topologies.idOf("T"));
            CHECK(reaction.top_type_to1() == topologies.idOf("TT"));
            CHECK(reaction.type_to1() == context.particleTypes().idOf("B"));
            CHECK(reaction.type_to2() == context.particleTypes().idOf("A"));
            CHECK_FALSE(reaction.is_topology_particle_reaction());
            CHECK(reaction.is_topology_topology_reaction());
            CHECK_FALSE(reaction.is_enzymatic());
            CHECK(reaction.is_fusion());
            CHECK(reaction.rate() == 10);
            CHECK(reaction.radius() == 20);
            CHECK(reaction.allow_self_connection());
            CHECK(reaction.mode() == readdy::model::top::reactions::STRMode::TT_FUSION_ALLOW_SELF);
        }
    }
    SECTION("Reactions") {
        context.particleTypes().add("A", 1.);
        context.particleTypes().add("B", 1.);
        context.particleTypes().add("C", 1.);

        SECTION("Decay") {
            context.reactions().add("mydecay:A->", 2.);
            const auto &r = context.reactions().order1ByName("mydecay");
            REQUIRE(r != nullptr);
            CHECK(r->type() == readdy::model::reactions::ReactionType::Decay);
            CHECK(r->nEducts() == 1);
            CHECK(r->nProducts() == 0);
            CHECK(r->educts()[0] == context.particleTypes().idOf("A"));
            CHECK(r->rate() == 2.);
        }
        SECTION("Conversion") {
            context.reactions().add("myconv: A -> B", 3.);
            const auto &r = context.reactions().order1ByName("myconv");
            REQUIRE(r != nullptr);
            CHECK(r->type() == readdy::model::reactions::ReactionType::Conversion);
            CHECK(r->nEducts() == 1);
            CHECK(r->nProducts() == 1);
            CHECK(r->educts()[0] == context.particleTypes().idOf("A"));
            CHECK(r->products()[0] == context.particleTypes().idOf("B"));
            CHECK(r->rate() == 3.);
        }
        SECTION("Fusion") {
            context.reactions().add("myfus: B +(1.2) B -> C [0.5, 0.5]", 4.);
            const auto &r = context.reactions().order2ByName("myfus");
            REQUIRE(r != nullptr);
            CHECK(r->type() == readdy::model::reactions::ReactionType::Fusion);
            CHECK(r->nEducts() == 2);
            CHECK(r->nProducts() == 1);
            CHECK(r->educts()[0] == context.particleTypes().idOf("B"));
            CHECK(r->educts()[1] == context.particleTypes().idOf("B"));
            CHECK(r->products()[0] == context.particleTypes().idOf("C"));
            CHECK(r->eductDistance() == 1.2);
            CHECK(r->weight1() == 0.5);
            CHECK(r->weight2() == 0.5);
            CHECK(r->rate() == 4.);
        }
        SECTION("Fission") {
            context.reactions().add("myfiss: B -> C +(3.0) B [0.1, 0.9]", 5.);
            const auto &r = context.reactions().order1ByName("myfiss");
            REQUIRE(r != nullptr);
            CHECK(r->type() == readdy::model::reactions::ReactionType::Fission);
            CHECK(r->nEducts() == 1);
            CHECK(r->nProducts() == 2);
            CHECK(r->educts()[0] == context.particleTypes().idOf("B"));
            CHECK(r->products()[0] == context.particleTypes().idOf("C"));
            CHECK(r->products()[1] == context.particleTypes().idOf("B"));
            CHECK(r->productDistance() == 3.0);
            CHECK(r->weight1() == 0.1);
            CHECK(r->weight2() == 0.9);
            CHECK(r->rate() == 5.);
        }
        SECTION("Enzymatic") {
            context.reactions().add("myenz:A +(1.5) C -> B + C", 6.);
            const auto &r = context.reactions().order2ByName("myenz");
            REQUIRE(r != nullptr);
            CHECK(r->type() == readdy::model::reactions::ReactionType::Enzymatic);
            CHECK(r->nEducts() == 2);
            CHECK(r->nProducts() == 2);
            CHECK(r->educts()[0] == context.particleTypes().idOf("A"));
            CHECK(r->educts()[1] == context.particleTypes().idOf("C"));
            CHECK(r->products()[0] == context.particleTypes().idOf("B"));
            CHECK(r->products()[1] == context.particleTypes().idOf("C"));
            CHECK(r->eductDistance() == 1.5);
            CHECK(r->rate() == 6.);
        }
        SECTION("Enzymatic predator-prey") {
            context.reactions().add("eat:B +(1) A -> A + A", 7.);
            const auto &r = context.reactions().order2ByName("eat");
            REQUIRE(r != nullptr);
            CHECK(r->type() == readdy::model::reactions::ReactionType::Enzymatic);
            CHECK(r->nEducts() == 2);
            CHECK(r->nProducts() == 2);
            CHECK(r->educts()[0] == context.particleTypes().idOf("B"));
            CHECK(r->educts()[1] == context.particleTypes().idOf("A"));
            CHECK(r->products()[0] == context.particleTypes().idOf("A"));
            CHECK(r->products()[1] == context.particleTypes().idOf("A"));
            CHECK(r->eductDistance() == 1.);
            CHECK(r->rate() == 7.);
        }

        SECTION("Invalid descriptors") {
            std::vector<std::string> inv = {
                    "myinvalid: + A -> B",
                    "noarrow: A B",
                    " : Noname ->",
                    "weights: A + A -> A [0.1, ]",
                    "blub: A (3)+ A -> B",
                    "eat: A +(1) B -> A + A"
            };
            for (const auto &i : inv) {
                REQUIRE_THROWS(context.reactions().add(i, 42.));
            }
            REQUIRE(context.reactions().nOrder1() == 0);
            REQUIRE(context.reactions().nOrder2() == 0);
        }

        SECTION("Existing reaction name") {
            context.reactions().add("bla: A->", 1.);
            REQUIRE_THROWS(context.reactions().add("bla: A->", 1.));
        }

        SECTION("Names and IDs") {
            context.reactions().add("foo: A->", 1.);
            context.reactions().add("bla: A+(1)A->A", 1.);
            const auto idFoo = context.reactions().idOf("foo");
            const auto idBla = context.reactions().idOf("bla");
            CHECK(idFoo >= 0);
            CHECK(idBla >= 0);
            CHECK(context.reactions().nameOf(idFoo) == "foo");
            CHECK(context.reactions().nameOf(idBla) == "bla");
        }

        SECTION("Fetching all reactions") {
            context.reactions().add("foo: A->", 1.);
            context.reactions().add("bla: A+(1)A->A", 1.);
            context.reactions().addConversion("conv1", "A", "A", 1.);
            context.reactions().addConversion("conv2", "A", "A", 1.);
            context.reactions().addFusion("fusion", "A","A", "A", 1., 1.);
            const auto &o1flat = context.reactions().order1Flat();
            const auto &o2flat = context.reactions().order2Flat();
            REQUIRE(o1flat.size() + o2flat.size() == 5);
        }
    }

    SECTION("Potentials") {
        context.particleTypes().add("a", 1.);
        context.particleTypes().add("b", 1.);
        auto noop = std::make_unique<readdy::testing::NOOPPotentialOrder2>(context.particleTypes()("a"),
                context.particleTypes()("b"));
        context.potentials().addUserDefined(noop.get());
        auto noop2 = std::make_unique<readdy::testing::NOOPPotentialOrder2>(
                context.particleTypes()("b"), context.particleTypes()("a"));
        context.potentials().addUserDefined(noop2.get());

        context.potentials().addScreenedElectrostatics("a", "b", 1., 1., 1., 1., 3, .1);
        {
            const auto &vector = context.potentials().potentialsOf("b", "a");
            REQUIRE(vector.size() == 3);
            CHECK(vector[0] == noop.get());
            CHECK(vector[1] == noop2.get());
        }
        {
            const auto &vector = context.potentials().potentialsOf("a", "b");
            REQUIRE(vector.size() == 3);
            CHECK(vector[0] == noop.get());
            CHECK(vector[1] == noop2.get());
        }
    }

    SECTION("Copyability") {
        Context ctx;
        ctx.particleTypes().add("A", 1.0);
        ctx.potentials().addHarmonicRepulsion("A", "A", 1.0, 1.0);
        ctx.reactions().addDecay("decay", "A", 1.0);
        ctx.topologyRegistry().addType("T");
        ctx.compartments().addSphere({{"A", "A"}}, "s",{0.,0.,0.}, 2., true);
        CHECK(ctx.particleTypes().idOf("A") >= 0);
        CHECK(ctx.particleTypes().nTypes() == 1);
        CHECK(ctx.potentials().potentialsOrder2().size() == 1);
        CHECK(ctx.potentials().potentialsOf("A", "A").size() == 1);
        CHECK(ctx.reactions().nOrder1() == 1);
        CHECK(ctx.reactions().idOf("decay") >= 0);
        CHECK(ctx.reactions().order1ByType("A").size() == 1);
        CHECK(ctx.topologyRegistry().idOf("T") >= 0);
        CHECK_FALSE(ctx.topologyRegistry().isSpatialReactionType("A"));
        CHECK(ctx.compartments().get().size() == 1);

        Context ctx2(ctx);
        CHECK(ctx2.particleTypes().idOf("A") >= 0);
        CHECK(ctx2.particleTypes().nTypes() == 1);
        CHECK(ctx2.potentials().potentialsOrder2().size() == 1);
        CHECK(ctx2.potentials().potentialsOf("A", "A").size() == 1);
        CHECK(ctx2.reactions().nOrder1() == 1);
        CHECK(ctx2.reactions().idOf("decay") >= 0);
        CHECK(ctx2.reactions().order1ByType("A").size() == 1);
        CHECK(ctx2.topologyRegistry().idOf("T") >= 0);
        CHECK_FALSE(ctx2.topologyRegistry().isSpatialReactionType("A"));
        CHECK(ctx2.compartments().get().size() == 1);
    }

    SECTION("Movability") {
        Context ctx;
        ctx.particleTypes().add("A", 1.0);
        ctx.potentials().addHarmonicRepulsion("A", "A", 1.0, 1.0);
        ctx.reactions().addDecay("decay", "A", 1.0);
        ctx.topologyRegistry().addType("T");
        ctx.compartments().addSphere({{"A", "A"}}, "s",{0.,0.,0.}, 2., true);
        CHECK(ctx.particleTypes().idOf("A") >= 0);
        CHECK(ctx.particleTypes().nTypes() == 1);
        CHECK(ctx.potentials().potentialsOrder2().size() == 1);
        CHECK(ctx.potentials().potentialsOf("A", "A").size() == 1);
        CHECK(ctx.reactions().nOrder1() == 1);
        CHECK(ctx.reactions().idOf("decay") >= 0);
        CHECK(ctx.reactions().order1ByType("A").size() == 1);
        CHECK(ctx.topologyRegistry().idOf("T") >= 0);
        CHECK_FALSE(ctx.topologyRegistry().isSpatialReactionType("A"));
        CHECK(ctx.compartments().get().size() == 1);

        Context ctx2(std::move(ctx));
        REQUIRE(&(ctx.particleTypes()) != &(ctx2.particleTypes())); // the reference has indeed been reset
        CHECK(ctx2.particleTypes().idOf("A") >= 0);
        CHECK(ctx2.particleTypes().nTypes() == 1);
        CHECK(ctx2.potentials().potentialsOrder2().size() == 1);
        CHECK(ctx2.potentials().potentialsOf("A", "A").size() == 1);
        CHECK(ctx2.reactions().nOrder1() == 1);
        CHECK(ctx2.reactions().idOf("decay") >= 0);
        CHECK(ctx2.reactions().order1ByType("A").size() == 1);
        CHECK(ctx2.topologyRegistry().idOf("T") >= 0);
        CHECK_FALSE(ctx2.topologyRegistry().isSpatialReactionType("A"));
        CHECK(ctx2.compartments().get().size() == 1);
    }

    SECTION("Adding reactions to the copy should not influence the original") {
        GIVEN("Context with one reaction for A, and a copy of that context") {
            Context ctx;
            ctx.particleTypes().add("A", 1.0);
            ctx.reactions().addDecay("decay", "A", 1.0);
            Context ctx2(ctx);
            WHEN("The copy gets another reaction for A") {
                ctx2.reactions().addFission("fiss", "A", "A", "A", 1.0, 1.0);
                THEN("The original context still only has one reaction for A") {
                    CHECK(ctx.reactions().order1ByType("A").size() == 1);
                }
            }
        }
    }

    SECTION("Modifying reactions of the copy should not influence the original") {
        Context ctx;
        ctx.particleTypes().add("A", 1.0);
        ctx.reactions().addDecay("decay", "A", 1.0);
        Context ctx2(ctx);
        WHEN("The copy's reaction's rate is modified") {
            REQUIRE(ctx2.reactions().order1Flat().size() == 1);
            ctx2.reactions().order1Flat().at(0)->rate() = 2.0;
            THEN("The original reaction still has a rate of 1.0") {
                REQUIRE(ctx.reactions().order1Flat().size() == 1);
                // todo what is the expected behavior here? deep copy or reference?

                // currently both contexts point to the same reaction due to shared_ptr
                // hence the following check fails
                // fixme CHECK(ctx.reactions().order1Flat().at(0)->rate() == 1.0);

                // reactions don't need polymorphism, the semantics are controlled by an enum anyway
                // i.e. the reaction registry can have vector<Reaction> (instead of vector<shared_ptr<Reaction>>)
                // then copy/move is trivial

                // potentials (and compartments) on the other hand need polymorphism due to virtual calculate...() functions
                // but registry could hold unique pointers and potentials implement deep copies
                // for deep copyability of context
            }
        }
    }

    SECTION("Test set kernel configuration") {
        Context ctx;
        WHEN("string is not valid json string (has missing braces)") {
            std::string invalid = R"({"MPI":{"dx":4.9,"dy":4.9,"dz":4.9,"haloThickness":1.0})";
            THEN("there must be an error") {
                REQUIRE_THROWS(ctx.setKernelConfiguration(invalid));
            }
        }
        WHEN("string is valid") {
            std::string valid = R"({"MPI":{"dx":4.9,"dy":5.9,"dz":6.9,"haloThickness":1.0}})";
            THEN("everything's OK and the appropriate values are set") {
                ctx.setKernelConfiguration(valid);
                auto& cfg = ctx.kernelConfiguration();
                REQUIRE(cfg.mpi.dx == Catch::Approx(4.9));
                REQUIRE(cfg.mpi.dy == Catch::Approx(5.9));
                REQUIRE(cfg.mpi.dz == Catch::Approx(6.9));
                REQUIRE(cfg.mpi.haloThickness == Catch::Approx(1.0));
            }
        }
    }
}
