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
 * Test different contexts w.r.t. boxsize and periodicity, perform setupBoxes() and see if that worked.
 * Then add a few particles and perform fillBoxes(). The result of this is the actual neighborlist, which is
 * checked as well.
 *
 * @file CPUTestNeighborList.cpp
 * @brief Test the neighborlist object of the CPU kernel.
 * @author chrisfroe
 * @author clonker
 * @date 23.08.16
 */

#include <cmath>

#include <gtest/gtest.h>
#include <readdy/api/SimulationScheme.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/cpu/nl/AdaptiveNeighborList.h>
#include <readdy/kernel/cpu/nl/CellDecompositionNeighborList.h>
#include <readdy/testing/NOOPPotential.h>


namespace cpu = readdy::kernel::cpu;
namespace m = readdy::model;

namespace {

using data_t = cpu::data::NLDataContainer;
using nl_t = cpu::cell_decomposition_neighbor_list;

struct TestNeighborList : ::testing::Test {

    std::unique_ptr<cpu::CPUKernel> kernel;
    readdy::particle_type_type typeIdA;

    std::unique_ptr<readdy::testing::NOOPPotentialOrder2> noop;

    TestNeighborList() : kernel(std::make_unique<cpu::CPUKernel>()) {
        auto &ctx = kernel->context();
        ctx.particle_types().add("A", 1.);
        auto a_id = ctx.particle_types()("A");
        readdy::scalar eductDistance = 1.2;
        ctx.reactions().addFusion("test", "A", "A", "A", 0., eductDistance);

        noop = std::make_unique<readdy::testing::NOOPPotentialOrder2>(a_id, a_id, 1.1, 0., 0.);
        ctx.potentials().addUserDefined(noop.get());
        typeIdA = ctx.particle_types().idOf("A");
        ctx.configure();
    }

};

class TestNeighborListImpl : public ::testing::TestWithParam<const char*> {};

auto isPairInList = [](readdy::kernel::cpu::nl::NeighborList *pairs, std::size_t idx1, std::size_t idx2) {
    bool foundOneDirection {false};
    bool foundOtherDirection {false};
    for(auto it : *pairs) {
        if(it.current_particle() == idx1) {
            foundOneDirection = std::find(it.begin(), it.end(), idx2) != it.end();
        } else {
            foundOtherDirection = std::find(it.begin(), it.end(), idx1) != it.end();
        }
        if(foundOneDirection && foundOtherDirection) {
            return true;
        }
    }
    return false;
};

auto isIdPairInList = [](nl_t *pairs, readdy::kernel::cpu::data::EntryDataContainer *data, std::size_t id1, std::size_t id2) {
    return isPairInList(pairs, data->getIndexForId(id1), data->getIndexForId(id2));
};

auto getNumberPairs = [](nl_t &pairs) {
    using val_t = decltype(*pairs.begin());
    return std::accumulate(pairs.begin(), pairs.end(), 0, [](int acc, val_t &x) {
        return acc + (x.end() - x.begin());
    });
};

TEST_F(TestNeighborList, ThreeBoxesNonPeriodic) {
    // maxcutoff is 1.2, system is 1.5 x 4 x 1.5, non-periodic, three cells
    auto &ctx = kernel->context();
    ctx.boxSize() = {{1.5, 4, 1.5}};
    ctx.periodicBoundaryConditions() = {{false, false, false}};

    readdy::util::thread::Config conf;

    nl_t list(ctx, conf);

    auto data = list.data();

    // Add three particles, two are in one outer box, the third on the other end and thus no neighbor
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, 1.8, 0, typeIdA)
    };

    data->addParticles(particles);
    list.set_up({});

    int sum = getNumberPairs(list);
    EXPECT_EQ(sum, 2);
    EXPECT_TRUE(isPairInList(&list, 0, 1));
    EXPECT_TRUE(isPairInList(&list, 1, 0));
}

TEST_F(TestNeighborList, OneDirection) {
    // maxcutoff is 1.2, system is 4.8 x 5 x 5.1
    auto &ctx = kernel->context();
    ctx.boxSize() = {{1.2, 1.1, 2.8}};
    ctx.periodicBoundaryConditions() = {{false, false, true}};
    ctx.potentials().addBox("A", .0, {-.4, -.4, -1.3}, {.4, .4, 1.3});
    ctx.configure(false);

    readdy::util::thread::Config conf;
    nl_t list(ctx, conf);
    // Add three particles, one of which is in the neighborhood of the other two
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, 0, -1.1, typeIdA), m::Particle(0, 0, .4, typeIdA), m::Particle(0, 0, 1.1, typeIdA)
    };
    std::vector<std::size_t> ids(particles.size());
    std::transform(particles.begin(), particles.end(), ids.begin(), [](const m::Particle& p) {return p.getId();});
    auto data = list.data();
    data->addParticles(particles);

    list.set_up({});

    int sum = getNumberPairs(list);
    EXPECT_LE(4, sum);
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(0), ids.at(2)));
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(2), ids.at(0)));
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(1), ids.at(2)));
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(2), ids.at(1)));
    EXPECT_FALSE(isIdPairInList(&list, data, ids.at(0), ids.at(1)));
    EXPECT_FALSE(isIdPairInList(&list, data, ids.at(1), ids.at(0)));
}

TEST_F(TestNeighborList, AllNeighborsInCutoffSphere) {
    // maxcutoff is 1.2, system is 4 x 4 x 4, all directions periodic
    auto &ctx = kernel->context();
    ctx.boxSize() = {{4, 4, 4}};
    ctx.periodicBoundaryConditions() = {{true, true, true}};
    readdy::util::thread::Config conf;
    nl_t list(ctx, conf);
    auto data = list.data();
    // Create a few particles. In this box setup, all particles are neighbors.
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0, 0, typeIdA), m::Particle(.3, 0, 0, typeIdA),
            m::Particle(0, .3, -.3, typeIdA), m::Particle(-.3, 0, .3, typeIdA), m::Particle(.3, -.3, 0, typeIdA)
    };

    data->addParticles(particles);
    list.set_up({});
    int sum = getNumberPairs(list);
    EXPECT_EQ(sum, 30);
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = i + 1; j < 6; ++j) {
            EXPECT_TRUE(isPairInList(&list, i, j)) << "Particles " << i << " and " << j << " were not neighbors.";
            EXPECT_TRUE(isPairInList(&list, j, i)) << "Particles " << j << " and " << i << " were not neighbors.";
        }
    }
}


TEST(TestAdaptiveNeighborList, CellContainerSanity) {
    using namespace readdy;
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    auto &context = kernel->context();
    context.boxSize() = {{10, 10, 10}};

    auto data = cpu::data::NLDataContainer(kernel->getCPUKernelStateModel().getParticleData());
    kernel::cpu::nl::CellContainer cellContainer{
            data,
            context,
            kernel->threadConfig()};

    ASSERT_EQ(cellContainer.size().x, 10);
    ASSERT_EQ(cellContainer.size().y, 10);
    ASSERT_EQ(cellContainer.size().z, 10);
    ASSERT_EQ(cellContainer.contiguous_index(), 0);

    cellContainer.subdivide(1);
    ASSERT_EQ(cellContainer.n_sub_cells()[0], 10);
    ASSERT_EQ(cellContainer.n_sub_cells()[1], 10);
    ASSERT_EQ(cellContainer.n_sub_cells()[2], 10);
    ASSERT_EQ(cellContainer.sub_cells().size(), 10 * 10 * 10);

    for (int i = 0; i < cellContainer.sub_cells().size(); ++i) {
        ASSERT_EQ(cellContainer.sub_cells().at(i).contiguous_index(), i);
    }
}

TEST(TestAdaptiveNeighborList, FirstLevelNeighborshipPeriodic) {
    using namespace readdy;
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    auto &context = kernel->context();
    context.periodicBoundaryConditions() = {{true, true, true}};
    context.boxSize() = {{10, 10, 10}};
    context.configure();

    auto data = data_t(kernel->getCPUKernelStateModel().getParticleData());
    kernel::cpu::nl::CellContainer cellContainer{
            data,
            context,
            kernel->threadConfig()};
    cellContainer.subdivide(1);
    cellContainer.refine_uniformly();
    cellContainer.setup_uniform_neighbors();

    {
        // corner case: lower left
        auto cell = cellContainer.leaf_cell_for_position({-4.9, -4.9, -4.9});
        ASSERT_TRUE(cell != nullptr);
        ASSERT_TRUE(cell->is_leaf());
        ASSERT_EQ(cell->level(), 2);
        ASSERT_TRUE(cell->offset() == readdy::Vec3(0, 0, 0));
        ASSERT_TRUE(cell->size() == readdy::Vec3(.5, .5, .5));
        ASSERT_EQ(cell->neighbors().size(), 5 * 5 * 5 - 1);
        std::set<readdy::Vec3> offsets;
        for (auto o : cell->neighbors()) {
            offsets.emplace(o->offset());
        }
        ASSERT_EQ(offsets.size(), 5 * 5 * 5 - 1);
    }
    {
        // cell in the middle
        auto cell = cellContainer.leaf_cell_for_position({0.1, 0.1, 0.1});
        ASSERT_TRUE(cell != nullptr);
        ASSERT_TRUE(cell->is_leaf());
        ASSERT_EQ(cell->level(), 2);
        ASSERT_TRUE(cell->offset() == readdy::Vec3(5, 5, 5));
        ASSERT_TRUE(cell->size() == readdy::Vec3(.5, .5, .5));
        const auto &cell_neighbors = cell->neighbors();
        ASSERT_EQ(cell_neighbors.size(), 5 * 5 * 5 - 1);
        auto shift = .5 * readdy::Vec3(context.boxSize());
        for (int i = -2; i <= 2; ++i) {
            for (int j = -2; j <= 2; ++j) {
                for (int k = -2; k <= 2; ++k) {
                    if (i == 0 && j == 0 && k == 0) continue;
                    auto pos = cell->offset() - shift +
                            readdy::Vec3(i * cell->size().x, j * cell->size().y, k * cell->size().z);
                    auto leaf = cellContainer.leaf_cell_for_position(pos);
                    ASSERT_TRUE(leaf != nullptr);
                    ASSERT_TRUE(std::find(cell_neighbors.begin(), cell_neighbors.end(), leaf) != cell_neighbors.end())
                                                << "Expect cell with lower left corner " << pos
                                                << " (resulting in cell with offset " << leaf->offset()
                                                << ") to be in the neighbors list";
                }
            }
        }
    }
}

TEST(TestAdaptiveNeighborList, FirstLevelNeighborshipNotPeriodic) {
    using namespace readdy;
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    auto &context = kernel->context();
    context.periodicBoundaryConditions() = {{false, false, false}};
    context.boxSize() = {{10, 10, 10}};

    auto data = data_t(kernel->getCPUKernelStateModel().getParticleData());
    kernel::cpu::nl::CellContainer cellContainer{
            data,
            context,
            kernel->threadConfig()};
    cellContainer.subdivide(1);
    cellContainer.refine_uniformly();
    cellContainer.setup_uniform_neighbors();

    {
        // corner case: lower left
        auto cell = cellContainer.leaf_cell_for_position({-4.9, -4.9, -4.9});
        ASSERT_TRUE(cell != nullptr);
        ASSERT_TRUE(cell->is_leaf());
        ASSERT_EQ(cell->level(), 2);
        ASSERT_TRUE(cell->offset() == readdy::Vec3(0, 0, 0));
        ASSERT_TRUE(cell->size() == readdy::Vec3(.5, .5, .5));
        const auto &cell_neighbors = cell->neighbors();
        ASSERT_EQ(cell_neighbors.size(), 3 * 3 * 3 - 1);
        {
            std::set<readdy::Vec3> offsets;
            for (auto o : cell_neighbors) {
                offsets.emplace(o->offset());
            }
            ASSERT_EQ(offsets.size(), 3 * 3 * 3 - 1);
        }
        {
            auto shift = .5 * Vec3(context.boxSize());
            for (int i = 0; i <= 2; ++i) {
                for (int j = 0; j <= 2; ++j) {
                    for (int k = 0; k <= 2; ++k) {
                        if (i == 0 && j == 0 && k == 0) continue;
                        auto pos = cell->offset() - shift +
                                   readdy::Vec3(i * cell->size().x, j * cell->size().y, k * cell->size().z);
                        auto leaf = cellContainer.leaf_cell_for_position(pos);
                        ASSERT_TRUE(leaf != nullptr);
                        ASSERT_TRUE(
                                std::find(cell_neighbors.begin(), cell_neighbors.end(), leaf) != cell_neighbors.end())
                                                    << "Expect cell with lower left corner " << pos
                                                    << " (resulting in cell with offset " << leaf->offset()
                                                    << ") to be in the neighbors list";
                    }
                }
            }
        }
    }
}

TEST(TestAdaptiveNeighborList, PartiallyPeriodicTube) {
    using namespace readdy;
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    auto &context = kernel->context();
    context.periodicBoundaryConditions() = {{false, true, true}};
    context.boxSize() = {{5, 1, 1.1}};
    context.configure();

    auto data = data_t(kernel->getCPUKernelStateModel().getParticleData());
    kernel::cpu::nl::CellContainer cellContainer{
            data,
            context,
            kernel->threadConfig()};
    cellContainer.subdivide(.5);
    cellContainer.refine_uniformly();
    cellContainer.setup_uniform_neighbors();

    EXPECT_EQ(cellContainer.sub_cells().size(), 10 * 2 * 2);
    {
        // upper right cell, check neighbors
        auto cell = cellContainer.leaf_cell_for_position({2.4, .4, .54});
        EXPECT_TRUE(cell != nullptr);
        EXPECT_EQ(cell->level(), 2);
        const auto &cell_neighbors = cell->neighbors();
        EXPECT_EQ(cell_neighbors.size(), 3 * 4 * 4 - 1);
        auto shift = .5 * Vec3(context.boxSize());
        for (int i = -2; i <= 0; ++i) {
            for (int j = -2; j <= 2; ++j) {
                for (int k = -2; k <= 2; ++k) {
                    if (i == 0 && j == 0 && k == 0) continue;
                    auto pos = context.applyPBCFun()(cell->offset() - shift +
                                                   readdy::Vec3(i * cell->size().x, j * cell->size().y,
                                                                       k * cell->size().z));
                    auto leaf = cellContainer.leaf_cell_for_position(pos);
                    ASSERT_TRUE(leaf != nullptr) << "should have a leaf for i=" << i << ", j=" << j << ", k=" << k
                                                 << ", pos=" << pos;
                    ASSERT_TRUE(
                            std::find(cell_neighbors.begin(), cell_neighbors.end(), leaf) != cell_neighbors.end())
                                                << "Expect cell with lower left corner " << pos
                                                << " (resulting in cell with offset " << leaf->offset()
                                                << ") to be in the neighbors list";
                }
            }
        }
    }
}

TEST(TestAdaptiveNeighborList, PositionToCell) {
    using namespace readdy;
    std::vector<Vec3> positions{{5,  5,  5},
                                       {0,  0,  0},
                                       {-5, -5, -5}};

    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    auto &context = kernel->context();
    context.periodicBoundaryConditions() = {{true, true, true}};
    context.boxSize() = {{10, 10, 10}};
    context.configure();

    auto data = data_t(kernel->getCPUKernelStateModel().getParticleData());
    kernel::cpu::nl::CellContainer cellContainer{
            data,
            context,
            kernel->threadConfig()};
    cellContainer.subdivide(1);
    cellContainer.refine_uniformly();
    cellContainer.setup_uniform_neighbors();

    const auto &pbc = context.applyPBCFun();
    auto shift = .5 * Vec3(context.boxSize());
    for (const auto &pos : positions) {
        auto proj_pos = pbc(pos);
        auto leaf = cellContainer.leaf_cell_for_position(proj_pos);
        ASSERT_TRUE(leaf != nullptr);
        ASSERT_TRUE(proj_pos >= leaf->offset() - shift);
        ASSERT_TRUE(proj_pos < leaf->offset() + leaf->size() - shift);
    }
}

TEST(TestAdaptiveNeighborList, SetUpNeighborList) {
    using namespace readdy;
    

    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();
    auto data = kernel->getCPUKernelStateModel().getParticleData();
    auto &context = kernel->context();
    const auto &pbc = context.applyPBCFun();

    context.periodicBoundaryConditions() = {{true, true, true}};
    context.boxSize() = {{10, 10, 10}};
    context.particle_types().add("A", 1.0);
    context.reactions().addFusion("test", "A", "A", "A", .1, .5);

    context.configure();

    auto n3 = readdy::model::rnd::normal3<readdy::scalar>;
    for (std::size_t i = 0; i < 500; ++i) {
        kernel->addParticle("A", pbc(n3(0, 10)));
    }

    kernel::cpu::nl::AdaptiveNeighborList neighbor_list{
            data,
            kernel->context(),
            kernel->threadConfig()
    };
    neighbor_list.set_up({});
    data = neighbor_list.data();

    for (std::size_t i = 0; i < data->size(); ++i) {
        auto cell = neighbor_list.cell_container().leaf_cell_for_position(pbc(data->pos(i)));
        ASSERT_TRUE(cell != nullptr) << "expected cell for position " << data->pos(i);
        const auto &cells_particles = cell->particles();
        ASSERT_TRUE(std::find(cells_particles.begin(), cells_particles.end(), i) != cells_particles.end());
    }
}

TEST(TestAdaptiveNeighborList, HilbertSort) {
    using namespace readdy;

    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();
    auto &context = kernel->context();
    context.boxSize() = {{1, 1, 1}};
    context.particle_types().add("A", 1.0);

    std::size_t i = 0;
    for (auto x = static_cast<scalar>(-.5); x < static_cast<scalar>(.5); x += static_cast<scalar>(.1)) {
        for (auto y = static_cast<scalar>(-.5); y < static_cast<scalar>(.5); y += static_cast<scalar>(.1)) {
            for (auto z = static_cast<scalar>(-.5); z < static_cast<scalar>(.5); z += static_cast<scalar>(.1)) {
                kernel->addParticle("A", {x, y, z});
                ++i;
            }
        }
    }
    auto &data = *kernel->getCPUKernelStateModel().getParticleData();
    ASSERT_EQ(data.size(), i);
    data.hilbertSort(.01);
    ASSERT_EQ(data.getNDeactivated(), 0);
    ASSERT_EQ(data.size(), i);

    /**
     * Can be plotted by
     *  import matplotlib as mpl
     *  from mpl_toolkits.mplot3d import Axes3D
     *  import numpy as np
     *  import matplotlib.pyplot as plt
     *  %matplotlib qt
     *
     *  hilbert = np.loadtxt("hilbert.txt", delimiter=",")
     *  X,Y,Z = hilbert[:,0], hilbert[:, 1], hilbert[:, 2]
     *
     *  fig = plt.figure()
     *  ax = fig.gca(projection='3d')
     *  ax.plot(X,Y,Z)
     *  ax.scatter(X,Y,Z)
     *  plt.show()
     *
     *  log::console()->set_pattern("%v");
     *  for(const auto& e : data) {
     *      log::debug("{}, {}, {}", e.position().x, e.position().y, e.position().z);
     *  }
     *
     */

    auto e0 = data.entry_at(0).id;
    auto e3 = data.entry_at(3).id;
    auto e15 = data.entry_at(15).id;

    data.removeEntry(0);
    data.removeEntry(3);
    data.removeEntry(15);

    ASSERT_EQ(data.getNDeactivated(), 3);
    data.hilbertSort(.01);
    data.blanks_moved_to_front();
    ASSERT_EQ(data.getNDeactivated(), 3);

    ASSERT_TRUE(data.entry_at(0).deactivated);
    ASSERT_TRUE(data.entry_at(1).deactivated);
    ASSERT_TRUE(data.entry_at(2).deactivated);

    ASSERT_TRUE(data.entry_at(0).id == e0 || data.entry_at(0).id == e3 || data.entry_at(0).id == e15);
    ASSERT_TRUE(data.entry_at(1).id == e0 || data.entry_at(1).id == e3 || data.entry_at(1).id == e15);
    ASSERT_TRUE(data.entry_at(2).id == e0 || data.entry_at(2).id == e3 || data.entry_at(2).id == e15);

    for (auto it = data.begin() + 3; it != data.end(); ++it) {
        ASSERT_FALSE(it->deactivated);
    }
}

TEST(TestAdaptiveNeighborList, VerletList) {
    using namespace readdy;

    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();
    auto &context = kernel->context();
    context.boxSize() = {{1, 1, 1}};
    context.particle_types().add("A", 1.0);

    for (int i = 0; i < 50; ++i) {
        kernel->addParticle("A", {static_cast<readdy::scalar>(model::rnd::uniform_real(-.5, .5)),
                                  static_cast<readdy::scalar>(model::rnd::uniform_real(-.5, .5)),
                                  static_cast<readdy::scalar>(model::rnd::uniform_real(-.5, .5))});
    }
    context.reactions().addFusion("test", "A", "A", "A", .1, .1);
    context.configure(false);
    const auto &d2 = context.distSquaredFun();
    auto data = kernel->getCPUKernelStateModel().getParticleData();
    kernel::cpu::nl::AdaptiveNeighborList neighbor_list{
            data,
            kernel->context(),
            kernel->threadConfig()
    };
    neighbor_list.set_up({});
    data = neighbor_list.data();

    std::size_t i = 0;
    for (const auto &entry_i : *data) {
        std::size_t j = 0;
        for (const auto &entry_j : *data) {
            if (!entry_i.deactivated && !entry_j.deactivated && i != j) {
                if (std::sqrt(d2(entry_i.pos, entry_j.pos)) < .1) {
                    const auto &neighbors = neighbor_list.nlData().neighbors().at(i);
                    ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), j) != neighbors.end())
                                                << i << " and " << j << " should be neighbors";
                }
            }
            ++j;
        }
        ++i;
    }
}

TEST(TestAdaptiveNeighborList, AdaptiveUpdating) {
    using namespace readdy;
    
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();
    auto &context = kernel->context();
    context.boxSize() = {{28, 28, 28}};
    context.particle_types().add("A", .1);
    context.particle_types().add("V", .1);
    context.periodicBoundaryConditions() = {{true, true, true}};
    context.kBT() = 1.0;

    auto cutoff = 1.5;
    for (int i = 0; i < 100; ++i) {
        kernel->addParticle("A", {static_cast<readdy::scalar>(model::rnd::uniform_real(-14., 14.)),
                                  static_cast<readdy::scalar>(model::rnd::uniform_real(-14., 14.)),
                                  static_cast<readdy::scalar>(model::rnd::uniform_real(-14., 14.))});
    }
    context.reactions().addFusion("test", "V", "V", "V", cutoff, cutoff);
    context.configure(false);
    const auto &d2 = context.distSquaredFun();
    auto data = kernel->getCPUKernelStateModel().getParticleData();
    kernel::cpu::nl::AdaptiveNeighborList neighbor_list{
            data,
            kernel->context(),
            kernel->threadConfig(),
            false
    };
    neighbor_list.skin() = 2.5;
    neighbor_list.set_up({});

    {
        std::size_t i = 0;
        for (const auto &entry_i : *data) {
            std::size_t j = 0;
            for (const auto &entry_j : *data) {
                if (!entry_i.deactivated && !entry_j.deactivated && i != j) {
                    if (d2(entry_i.pos, entry_j.pos) < cutoff * cutoff) {
                        const auto &neighbors = neighbor_list.nlData().neighbors().at(i);
                        ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), j) != neighbors.end())
                                                    << i << " and " << j << " should be neighbors";
                    }
                }
                ++j;
            }
            ++i;
        }
    }

    auto integrator = kernel->createAction<readdy::model::actions::EulerBDIntegrator>(.01);

    for (int t = 0; t < 100; ++t) {
        integrator->perform();
        neighbor_list.update({});
        {
            std::size_t i = 0;
            for (const auto &entry_i : *data) {
                std::size_t j = 0;
                for (const auto &entry_j : *data) {
                    if ((!entry_i.deactivated) && (!entry_j.deactivated) && i != j) {
                        if (d2(entry_i.pos, entry_j.pos) < cutoff * cutoff) {
                            const auto &neighbors = neighbor_list.nlData().neighbors().at(i);
                            bool neighbors_of_each_other = std::find(neighbors.begin(), neighbors.end(), j) != neighbors.end();
                            EXPECT_TRUE(neighbors_of_each_other);
                        }
                    }
                    ++j;
                }
                ++i;
            }
        }
    }


}


TEST_P(TestNeighborListImpl, DiffusionAndReaction) {
    using namespace readdy;
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    // A is absorbed and created by F, while the number of F stays constant, this test spans multiple timesteps
    kernel->context().particle_types().add("A", 0.05);
    kernel->context().particle_types().add("F", 0.0);
    kernel->context().particle_types().add("V", 0.0);
    kernel->context().periodicBoundaryConditions() = {{true, true, true}};
    kernel->context().boxSize() = {{100, 10, 10}};
    readdy::conf::Configuration conf {};
    conf.cpu.neighborList.type = GetParam();
    kernel->context().kernelConfiguration() = conf;

    const auto weightF = static_cast<readdy::scalar>(0.);
    const auto weightA = static_cast<readdy::scalar>(1.);
    kernel->context().reactions().addFusion("F+A->F", "A", "F", "F", .1, 2.0, weightF, weightA);
    //kernel->registerReaction<readdy::model::reactions::Fusion>("F+A->F2", "V", "V", "V", .1, 2.0, weightF, weightA);

    auto n3 = readdy::model::rnd::normal3<readdy::scalar>;
    // 120 F particles
    for (std::size_t i = 0; i < 100; ++i) {
        kernel->addParticle("F", n3(.0, 1.));
        kernel->addParticle("A", n3(0., 1.));
    }

    auto obs = kernel->createObservable<readdy::model::observables::NParticles>(1, std::vector<std::string>({"F", "A"}));
    obs->setCallback(
            [&](const readdy::model::observables::NParticles::result_type &result) {
                EXPECT_EQ(result[0], 100);
            }
    );
    auto connection = kernel->connectObservable(obs.get());

    {
        readdy::util::PerformanceNode pn("", false);
        auto conf = readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>(kernel.get(), pn);
        conf.withReactionScheduler<readdy::model::actions::reactions::Gillespie>();
        conf.withSkinSize(.1);
        conf.configureAndRun(100, .01);
    }
}

TEST_P(TestNeighborListImpl, Diffusion) {
    using namespace readdy;
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    // A is absorbed and created by F, while the number of F stays constant, this test spans multiple timesteps
    auto& context = kernel->context();
    context.particle_types().add("A", 0.05);
    context.particle_types().add("F", 0.0);
    context.particle_types().add("V", 0.0);

    scalar cutoff = 2.0;

    // just to have a cutoff
    context.reactions().addFusion("test", "V", "V", "V", .1, cutoff);
    context.periodicBoundaryConditions() = {{true, true, true}};
    context.boxSize() = {{100, 10, 10}};

    {
        readdy::conf::Configuration conf{};
        conf.cpu.neighborList.type = GetParam();
        kernel->context().kernelConfiguration() = conf;
    }

    auto n3 = readdy::model::rnd::normal3<readdy::scalar>;
    // 120 F particles
    for (std::size_t i = 0; i < 50; ++i) {
        kernel->addParticle("F", n3(0., 1.));
        kernel->addParticle("A", n3(0., 1.));
    }
    auto obs = kernel->createObservable<readdy::model::observables::NParticles>(1);
    obs->setCallback(
            [&](const readdy::model::observables::NParticles::result_type &) {
                const auto &d2 = context.distSquaredFun();
                const auto neighbor_list = kernel->getCPUKernelStateModel().getNeighborList();

                for(auto it = neighbor_list->begin(); it != neighbor_list->end(); ++it) {
                    const auto &entry = neighbor_list->data()->entry_at(it->current_particle());
                    ASSERT_FALSE(entry.deactivated);

                    std::vector<std::size_t> neighbors;
                    neighbors.reserve(it->n_neighbors());
                    for (unsigned long neighbor : *it) {
                        const auto &neighborEntry = neighbor_list->data()->entry_at(neighbor);
                        ASSERT_FALSE(neighborEntry.deactivated);
                        neighbors.push_back(neighbor);
                        ASSERT_LE(d2(entry.pos, neighborEntry.pos), (cutoff+.1) * (cutoff+.1));
                    }

                    std::size_t pidx = 0;
                    for(const auto &e : *neighbor_list->data()) {
                        ASSERT_FALSE(e.deactivated);
                        if (pidx != it->current_particle() && d2(entry.pos, e.pos) < cutoff * cutoff) {
                            ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                        }
                        ++pidx;
                    }
                }
            }
    );
    auto connection = kernel->connectObservable(obs.get());
    {
        readdy::util::PerformanceNode pn("", false);
        auto sc = readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>(kernel.get(), pn);
        sc.withReactionScheduler<readdy::model::actions::reactions::Gillespie>();
        sc.withSkinSize(.1);
        sc.configureAndRun(100, .01);
    }
}

}

INSTANTIATE_TEST_CASE_P(TestNeighborListImplImplementation, TestNeighborListImpl,
                        ::testing::Values("CompactCLL", "ContiguousCLL", "DynamicCLL", "CellDecomposition", "Adaptive"));
