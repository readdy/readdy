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
 * @file TestNeighborList2.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 24.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <gtest/gtest.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/cpu/nl/NeighborList.h>

namespace {

TEST(TestNeighborList2, CellContainerSanity) {
    using namespace readdy;
    log::console()->set_level(spdlog::level::debug);
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    auto &context = kernel->getKernelContext();
    context.setBoxSize(10, 10, 10);

    kernel::cpu::nl::CellContainer cellContainer{
            *kernel->getCPUKernelStateModel().getParticleData(),
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

TEST(TestNeighborList2, FirstLevelNeighborshipPeriodic) {
    using namespace readdy;
    log::console()->set_level(spdlog::level::debug);
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    auto &context = kernel->getKernelContext();
    context.setPeriodicBoundary(true, true, true);
    context.setBoxSize(10, 10, 10);

    kernel::cpu::nl::CellContainer cellContainer{
            *kernel->getCPUKernelStateModel().getParticleData(),
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
        ASSERT_TRUE(cell->offset() == readdy::model::Vec3(0, 0, 0));
        ASSERT_TRUE(cell->size() == readdy::model::Vec3(.5, .5, .5));
        ASSERT_EQ(cell->neighbors().size(), 5 * 5 * 5 - 1);
        std::set<readdy::model::Vec3> offsets;
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
        ASSERT_TRUE(cell->offset() == readdy::model::Vec3(5, 5, 5));
        ASSERT_TRUE(cell->size() == readdy::model::Vec3(.5, .5, .5));
        const auto &cell_neighbors = cell->neighbors();
        ASSERT_EQ(cell_neighbors.size(), 5 * 5 * 5 - 1);
        auto shift = .5 * model::Vec3(context.getBoxSize());
        for (int i = -2; i <= 2; ++i) {
            for (int j = -2; j <= 2; ++j) {
                for (int k = -2; k <= 2; ++k) {
                    if (i == 0 && j == 0 && k == 0) continue;
                    auto pos = cell->offset() - shift +
                               readdy::model::Vec3(i * cell->size().x, j * cell->size().y, k * cell->size().z);
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

TEST(TestNeighborList2, FirstLevelNeighborshipNotPeriodic) {
    using namespace readdy;
    log::console()->set_level(spdlog::level::debug);
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    auto &context = kernel->getKernelContext();
    context.setPeriodicBoundary(false, false, false);
    context.setBoxSize(10, 10, 10);

    kernel::cpu::nl::CellContainer cellContainer{
            *kernel->getCPUKernelStateModel().getParticleData(),
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
        ASSERT_TRUE(cell->offset() == readdy::model::Vec3(0, 0, 0));
        ASSERT_TRUE(cell->size() == readdy::model::Vec3(.5, .5, .5));
        const auto &cell_neighbors = cell->neighbors();
        ASSERT_EQ(cell_neighbors.size(), 3 * 3 * 3 - 1);
        {
            std::set<readdy::model::Vec3> offsets;
            for (auto o : cell_neighbors) {
                offsets.emplace(o->offset());
            }
            ASSERT_EQ(offsets.size(), 3 * 3 * 3 - 1);
        }
        {
            auto shift = .5 * model::Vec3(context.getBoxSize());
            for (int i = 0; i <= 2; ++i) {
                for (int j = 0; j <= 2; ++j) {
                    for (int k = 0; k <= 2; ++k) {
                        if (i == 0 && j == 0 && k == 0) continue;
                        auto pos = cell->offset() - shift +
                                   readdy::model::Vec3(i * cell->size().x, j * cell->size().y, k * cell->size().z);
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

TEST(TestNeighborList2, PartiallyPeriodicTube) {
    using namespace readdy;
    log::console()->set_level(spdlog::level::debug);
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    auto &context = kernel->getKernelContext();
    context.setPeriodicBoundary(false, true, true);
    context.setBoxSize(5, 1, 1.1);
    kernel::cpu::nl::CellContainer cellContainer{
            *kernel->getCPUKernelStateModel().getParticleData(),
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
        const auto &cell_neighbors = cell->neighbors();
        EXPECT_EQ(cell_neighbors.size(), 3 * 5 * 5 - 1);
        auto shift = .5 * model::Vec3(context.getBoxSize());
        for (int i = -2; i <= 0; ++i) {
            for (int j = -2; j <= 2; ++j) {
                for (int k = -2; k <= 2; ++k) {
                    if (i == 0 && j == 0 && k == 0) continue;
                    auto pos = context.getPBCFun()(cell->offset() - shift +
                                                   readdy::model::Vec3(i * cell->size().x, j * cell->size().y,
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

TEST(TestNeighborList2, PositionToCell) {
    using namespace readdy;
    std::vector<model::Vec3> positions{{5,  5,  5},
                                       {0,  0,  0},
                                       {-5, -5, -5}};

    log::console()->set_level(spdlog::level::debug);
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    auto &context = kernel->getKernelContext();
    context.setPeriodicBoundary(true, true, true);
    context.setBoxSize(10, 10, 10);
    kernel::cpu::nl::CellContainer cellContainer{
            *kernel->getCPUKernelStateModel().getParticleData(),
            context,
            kernel->threadConfig()};
    cellContainer.subdivide(1);
    cellContainer.refine_uniformly();
    cellContainer.setup_uniform_neighbors();

    const auto &pbc = context.getPBCFun();
    auto shift = .5 * model::Vec3(context.getBoxSize());
    for (const auto &pos : positions) {
        auto proj_pos = pbc(pos);
        auto leaf = cellContainer.leaf_cell_for_position(proj_pos);
        ASSERT_TRUE(leaf != nullptr);
        ASSERT_TRUE(proj_pos >= leaf->offset() - shift);
        ASSERT_TRUE(proj_pos < leaf->offset() + leaf->size() - shift);
    }
}

TEST(TestNeighborList2, SetUpNeighborList) {
    using namespace readdy;
    log::console()->set_level(spdlog::level::debug);

    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();
    const auto& data = *kernel->getCPUKernelStateModel().getParticleData();
    auto &context = kernel->getKernelContext();
    const auto& pbc = context.getPBCFun();

    context.setPeriodicBoundary(true, true, true);
    context.setBoxSize(10, 10, 10);
    context.particle_types().add("A", 1.0, 1.0);
    kernel->registerReaction<readdy::model::reactions::Fusion>("test", "A", "A", "A", .1, .5);

    context.configure();

    auto n3 = readdy::model::rnd::normal3<>;
    for (std::size_t i = 0; i < 500; ++i) {
        kernel->addParticle("A", pbc(n3(0, 10)));
    }

    kernel::cpu::nl::NeighborList neighbor_list {
            data,
            kernel->getKernelContext(),
            kernel->threadConfig()
    };
    neighbor_list.set_up();

    for(std::size_t i = 0; i < data.size(); ++i) {
        auto cell = neighbor_list.cell_container().leaf_cell_for_position(pbc(data.pos(i)));
        ASSERT_TRUE(cell != nullptr) << "expected cell for position " << data.pos(i);
        const auto& cells_particles = cell->particles().get();
        ASSERT_TRUE(std::find(cells_particles.begin(), cells_particles.end(), i) != cells_particles.end());
    }
}

TEST(TestNeighborList2, HilbertSort) {
    using namespace readdy;
    log::console()->set_level(spdlog::level::debug);

    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();
    auto &context = kernel->getKernelContext();
    context.setBoxSize(1, 1, 1);
    context.particle_types().add("A", 1.0, 1.0);

    for(scalar x = -.5; x < .5; x += .1) {
        for(scalar y = -.5; y < .5; y += .1) {
            for(scalar z = -.5; z < .5; z += .1) {
                kernel->addParticle("A", {x,y,z});
            }
        }
    }
    auto& data = *kernel->getCPUKernelStateModel().getParticleData();
    data.hilbert_sort(.01);
    ASSERT_EQ(data.getNDeactivated(), 0);
    ASSERT_EQ(data.size(), 10*10*10);

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
    auto e3= data.entry_at(3).id;
    auto e15 = data.entry_at(15).id;

    data.removeEntry(0);
    data.removeEntry(3);
    data.removeEntry(15);

    ASSERT_EQ(data.getNDeactivated(), 3);
    data.hilbert_sort(.01);
    data.blanks_moved_to_front();
    ASSERT_EQ(data.getNDeactivated(), 3);

    ASSERT_TRUE(data.entry_at(0).is_deactivated());
    ASSERT_TRUE(data.entry_at(1).is_deactivated());
    ASSERT_TRUE(data.entry_at(2).is_deactivated());

    ASSERT_TRUE(data.entry_at(0).id == e0 || data.entry_at(0).id == e3 || data.entry_at(0).id == e15);
    ASSERT_TRUE(data.entry_at(1).id == e0 || data.entry_at(1).id == e3 || data.entry_at(1).id == e15);
    ASSERT_TRUE(data.entry_at(2).id == e0 || data.entry_at(2).id == e3 || data.entry_at(2).id == e15);

    for(auto it = data.begin()+3; it != data.end(); ++it) {
        ASSERT_FALSE(it->is_deactivated());
    }
}

}