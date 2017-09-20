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
 * Testfile containing tests for the KernelContext class.
 *
 * @file TestKernelContext.cpp
 * @brief Testfile for KernelContext.
 * @author clonker
 * @date 19.04.16
 */


#include <readdy/common/Utils.h>
#include <readdy/common/make_unique.h>
#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/testing/NOOPPotential.h>
#include <readdy/model/potentials/PotentialsOrder1.h>

namespace m = readdy::model;

namespace {

class TestKernelContext : public ::testing::Test {
protected:
    TestKernelContext() {}
};

class TestKernelContextWithKernels : public KernelTest {

};

TEST_F(TestKernelContext, SetGetKBT) {
    m::KernelContext ctx;
    ctx.kBT() = 42;
    EXPECT_EQ(42, ctx.kBT());
}

TEST_F(TestKernelContext, PeriodicBoundary) {
    m::KernelContext ctx;
    ctx.periodicBoundaryConditions() = {{true, false, true}};
    auto boundary = ctx.periodicBoundaryConditions();
    EXPECT_TRUE(boundary[0]);
    EXPECT_FALSE(boundary[1]);
    EXPECT_TRUE(boundary[2]);
}

TEST_F(TestKernelContext, BoxSize) {
    m::KernelContext ctx;
    ctx.boxSize() = {{10, 11, 12}};
    auto box_size = ctx.boxSize();
    EXPECT_EQ(box_size[0], 10);
    EXPECT_EQ(box_size[1], 11);
    EXPECT_EQ(box_size[2], 12);
}

TEST_F(TestKernelContext, Copyable) {
    m::KernelContext context;
    m::KernelContext copy(context);
}

TEST_F(TestKernelContext, PotentialOrder2Map) {
    m::KernelContext ctx;
    ctx.particle_types().add("a", 1., 1.);
    ctx.particle_types().add("b", 1., 1.);
    auto noop = std::make_unique<readdy::testing::NOOPPotentialOrder2>(ctx.particle_types()("a"), ctx.particle_types()("b"));
    ctx.potentials().addUserDefined(noop.get());
    auto noop2 = std::make_unique<readdy::testing::NOOPPotentialOrder2>(ctx.particle_types()("b"),ctx.particle_types()( "a"));
    ctx.potentials().addUserDefined(noop2.get());
    ctx.configure();
    auto vector = ctx.potentials().potentials_of("b", "a");
    EXPECT_EQ(vector.size(), 2);
}

TEST_P(TestKernelContextWithKernels, PotentialOrder1Map) {
    using vec_t = readdy::Vec3;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");

    namespace rmp = readdy::model::potentials;

    auto &ctx = kernel->getKernelContext();

    ctx.particle_types().add("A", 1.0, 1.0);
    ctx.particle_types().add("B", 3.0, 2.0);
    ctx.particle_types().add("C", 4.0, 3.0);
    ctx.particle_types().add("D", 2.0, 4.0);

    std::vector<short> idsToRemove;
    short uuid2_1, uuid2_2;
    {
        auto id1 = kernel->getKernelContext().potentials().addCube("A", 0, vec_t{0, 0, 0}, vec_t{0, 0, 0});
        auto id2 = kernel->getKernelContext().potentials().addCube("C", 0, vec_t{0, 0, 0}, vec_t{0, 0, 0});
        auto id3 = kernel->getKernelContext().potentials().addCube("D", 0, vec_t{0, 0, 0}, vec_t{0, 0, 0});
        auto id4 = kernel->getKernelContext().potentials().addCube("C", 0, vec_t{0, 0, 0}, vec_t{0, 0, 0});

        idsToRemove.push_back(id1);
        idsToRemove.push_back(id2);
        idsToRemove.push_back(id3);
        idsToRemove.push_back(id4);

        kernel->getKernelContext().potentials().addCube("B", 0, vec_t{0, 0, 0}, vec_t{0, 0, 0});

        uuid2_1 = kernel->getKernelContext().potentials().addHarmonicRepulsion("A", "C", 0);
        uuid2_2 = kernel->getKernelContext().potentials().addHarmonicRepulsion("B", "C", 0);
        ctx.configure();
    }
    // test that order 1 potentials are set up correctly
    {
        {
            const auto &pot1_A = ctx.potentials().potentials_of("A");
            EXPECT_EQ(pot1_A.size(), 1);
            EXPECT_EQ(dynamic_cast<rmp::Cube *>(pot1_A[0])->getParticleRadius(), 1.0);
        }
        {
            const auto &pot1_B = ctx.potentials().potentials_of("B");
            EXPECT_EQ(pot1_B.size(), 1);
            EXPECT_EQ(dynamic_cast<rmp::Cube *>(pot1_B[0])->getParticleRadius(), 2.0);
        }
        {
            const auto &pot1_C = ctx.potentials().potentials_of("C");
            EXPECT_EQ(pot1_C.size(), 2);
            for (auto &&ptr : pot1_C) {
                EXPECT_EQ(dynamic_cast<rmp::Cube *>(ptr)->getParticleRadius(), 3.0);
            }
        }
        {
            const auto &pot1_D = ctx.potentials().potentials_of("D");
            EXPECT_EQ(pot1_D.size(), 1);
            EXPECT_EQ(dynamic_cast<rmp::Cube *>(pot1_D[0])->getParticleRadius(), 4.0);
        }
    }
    // test that order 2 potentials are set up correctly
    {
        EXPECT_EQ(ctx.potentials().potentials_of("A", "A").size(), 0);
        EXPECT_EQ(ctx.potentials().potentials_of("A", "B").size(), 0);
        EXPECT_EQ(ctx.potentials().potentials_of("A", "D").size(), 0);
        EXPECT_EQ(ctx.potentials().potentials_of("B", "B").size(), 0);
        EXPECT_EQ(ctx.potentials().potentials_of("B", "D").size(), 0);

        EXPECT_EQ(ctx.potentials().potentials_of("A", "C").size(), 1);
        EXPECT_EQ(ctx.potentials().potentials_of("B", "C").size(), 1);
        {
            const auto &pot2_AC = ctx.potentials().potentials_of("A", "C");
            EXPECT_EQ(pot2_AC[0]->getId(), uuid2_1);
            EXPECT_EQ(dynamic_cast<rmp::HarmonicRepulsion *>(pot2_AC[0])->getSumOfParticleRadii(), 1 + 3);
        }
        {
            const auto &pot2_BC = ctx.potentials().potentials_of("B", "C");
            EXPECT_EQ(pot2_BC[0]->getId(), uuid2_2);
            EXPECT_EQ(dynamic_cast<rmp::HarmonicRepulsion *>(pot2_BC[0])->getSumOfParticleRadii(), 2 + 3);
        }
    }

    // now remove
    std::for_each(idsToRemove.begin(), idsToRemove.end(), [&](const short id) { ctx.potentials().remove(id); });
    {
        // only one potential for particle type B has a different uuid
        ctx.configure();
        EXPECT_EQ(ctx.potentials().potentials_of("A").size(), 0);
        EXPECT_EQ(ctx.potentials().potentials_of("B").size(), 1);
        EXPECT_EQ(ctx.potentials().potentials_of("C").size(), 0);
        EXPECT_EQ(ctx.potentials().potentials_of("D").size(), 0);

        // remove 2nd potential
        ctx.potentials().remove(uuid2_2);
        ctx.configure();
        EXPECT_EQ(ctx.potentials().potentials_of("A", "C").size(), 1);
        EXPECT_EQ(ctx.potentials().potentials_of("B", "C").size(), 0);
    }

}

TEST_F(TestKernelContext, ReactionDescriptorAddReactions) {
    m::KernelContext ctx;
    ctx.particle_types().add("A", 1., 1.);
    ctx.particle_types().add("B", 1., 1.);
    ctx.particle_types().add("C", 1., 1.);
    {
        auto decay = "mydecay:A->";
        ctx.reactions().add(decay, 2.);
        const auto &r = ctx.reactions().order1_by_name("mydecay");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->getType(), m::reactions::ReactionType::Decay);
        EXPECT_EQ(r->getNEducts(), 1);
        EXPECT_EQ(r->getNProducts(), 0);
        EXPECT_EQ(r->getEducts()[0], ctx.particle_types().id_of("A"));
        EXPECT_EQ(r->getRate(), 2.);
    }
    {
        auto conversion = "myconv: A -> B";
        ctx.reactions().add(conversion, 3.);
        const auto &r = ctx.reactions().order1_by_name("myconv");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->getType(), m::reactions::ReactionType::Conversion);
        EXPECT_EQ(r->getNEducts(), 1);
        EXPECT_EQ(r->getNProducts(), 1);
        EXPECT_EQ(r->getEducts()[0], ctx.particle_types().id_of("A"));
        EXPECT_EQ(r->getProducts()[0], ctx.particle_types().id_of("B"));
        EXPECT_EQ(r->getRate(), 3.);
    }
    {
        auto fusion = "myfus: B +(1.2) B -> C [0.5, 0.5]";
        ctx.reactions().add(fusion, 4.);
        const auto &r = ctx.reactions().order2_by_name("myfus");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->getType(), m::reactions::ReactionType::Fusion);
        EXPECT_EQ(r->getNEducts(), 2);
        EXPECT_EQ(r->getNProducts(), 1);
        EXPECT_EQ(r->getEducts()[0], ctx.particle_types().id_of("B"));
        EXPECT_EQ(r->getEducts()[1], ctx.particle_types().id_of("B"));
        EXPECT_EQ(r->getProducts()[0], ctx.particle_types().id_of("C"));
        EXPECT_EQ(r->getEductDistance(), 1.2);
        EXPECT_EQ(r->getWeight1(), 0.5);
        EXPECT_EQ(r->getWeight2(), 0.5);
        EXPECT_EQ(r->getRate(), 4.);
    }
    {
        auto fission = "myfiss: B -> C +(3.0) B [0.1, 0.9]";
        ctx.reactions().add(fission, 5.);
        const auto &r = ctx.reactions().order1_by_name("myfiss");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->getType(), m::reactions::ReactionType::Fission);
        EXPECT_EQ(r->getNEducts(), 1);
        EXPECT_EQ(r->getNProducts(), 2);
        EXPECT_EQ(r->getEducts()[0], ctx.particle_types().id_of("B"));
        EXPECT_EQ(r->getProducts()[0], ctx.particle_types().id_of("C"));
        EXPECT_EQ(r->getProducts()[1], ctx.particle_types().id_of("B"));
        EXPECT_EQ(r->getProductDistance(), 3.0);
        EXPECT_EQ(r->getWeight1(), 0.1);
        EXPECT_EQ(r->getWeight2(), 0.9);
        EXPECT_EQ(r->getRate(), 5.);
    }
    {
        auto enzymatic = "myenz:A +(1.5) C -> B + C";
        ctx.reactions().add(enzymatic, 6.);
        const auto &r = ctx.reactions().order2_by_name("myenz");
        ASSERT_NE(r, nullptr);
        EXPECT_EQ(r->getType(), m::reactions::ReactionType::Enzymatic);
        EXPECT_EQ(r->getNEducts(), 2);
        EXPECT_EQ(r->getNProducts(), 2);
        EXPECT_EQ(r->getEducts()[0], ctx.particle_types().id_of("A"));
        EXPECT_EQ(r->getEducts()[1], ctx.particle_types().id_of("C"));
        EXPECT_EQ(r->getProducts()[0], ctx.particle_types().id_of("B"));
        EXPECT_EQ(r->getProducts()[1], ctx.particle_types().id_of("C"));
        EXPECT_EQ(r->getEductDistance(), 1.5);
        EXPECT_EQ(r->getRate(), 6.);
    }
}

TEST_F(TestKernelContext, ReactionDescriptorInvalidInputs) {
    m::KernelContext ctx;
    ctx.particle_types().add("A", 1., 1.);
    ctx.particle_types().add("B", 1., 1.);
    std::vector<std::string> inv = {"myinvalid: + A -> B", "noarrow: A B", " : Noname ->", "weights: A + A -> A [0.1, ]"};
    for (const auto &i : inv) {
        auto addInv = [&ctx, &i]() {
            ctx.reactions().add(i, 42.);
        };
        EXPECT_ANY_THROW(addInv());
    }
    EXPECT_EQ(ctx.reactions().n_order1(), 0);
    EXPECT_EQ(ctx.reactions().n_order2(), 0);
}

INSTANTIATE_TEST_CASE_P(TestKernelContext, TestKernelContextWithKernels,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}