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
 * @file TestHilbertCurveIndexing.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 10.11.16
 */

#include "gtest/gtest.h"

#include <hilbert.h>
#include <readdy/plugin/KernelProvider.h>

namespace {


struct Position {
    std::array<readdy::scalar, 3> data;

    Position(const std::array<readdy::scalar, 3> &data) : data(data) {}
    Position() : Position({0,0,0}) {}

    Position(Position&&) = default;
    Position& operator=(Position&&) = default;
    Position(const Position&) = delete;
    Position& operator=(const Position&) = delete;

};


/**
 * Plan:
 *  0. change data structure to vector of vectors
 *  1. Sort boxes
 *  2. assign particles to boxes, group
 *  3. after some time, periodically, reorganize the data incrementally
 *     (depending on (#(particles ungrouped)/#(particles total))?)
 */


TEST(TestCurveIndexing, Cells) {

    auto uniform_real = &readdy::model::rnd::uniform_real<>;

    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU_Legacy");
    kernel->context().boxSize() = {{5, 5, 5}};
    const auto &simBoxSize = kernel->context().boxSize();

    // 2x2x2 = 8 boxes
    const std::array<unsigned int, 3> nCells{{2, 2, 2}};
    const std::array<readdy::scalar, 3> cellWidth{{1, 1, 1}};


    const std::size_t data_size = 100;
    std::vector<int> hilbert_indices;
    hilbert_indices.resize(data_size);
    std::vector<std::size_t> indices(data_size);
    std::iota(indices.begin(), indices.end(), 0);

    {
        // assert indices are strictly increasing
        std::size_t old_idx = 0;
        for (const auto idx : indices) {
            ASSERT_TRUE(old_idx == 0 || idx == old_idx + 1);
            old_idx = idx;
        }
    }

    std::vector<Position> positions(data_size);
    for (auto &pos : positions) {
        pos.data[0] = uniform_real(0, 5);
        pos.data[1] = uniform_real(0, 15);
        pos.data[2] = uniform_real(0, 25);
    }

    {
        // fill up hilbert indices
        using cell_index = unsigned int;
        auto hilbert_it = hilbert_indices.begin();
        for (int _i = 0; _i < data_size; ++_i) {
            const auto &pos = positions[_i];
            const cell_index i = static_cast<const cell_index>(floor(pos.data[0] / .001));
            const cell_index j = static_cast<const cell_index>(floor(pos.data[1] / .001));
            const cell_index k = static_cast<const cell_index>(floor(pos.data[2] / .001));
            bitmask_t coords[3]{i, j, k};
            *hilbert_it = static_cast<unsigned int>(hilbert_c2i(3, CHAR_BIT, coords));
            ++hilbert_it;
        }
    }

    const auto n_threads = 8;
    const auto grainSize = 100 / n_threads;
    {
        {
            // sort
            std::sort(indices.begin(), indices.end(),
                      [&hilbert_indices](std::size_t i, std::size_t j) {
                          return hilbert_indices[i] < hilbert_indices[j];
                      });
        }

    }

    {
        int last_idx = -1;
        for (const auto sort_idx : indices) {
            ASSERT_TRUE(last_idx <= hilbert_indices[sort_idx]);
            last_idx = hilbert_indices[sort_idx];
        }
    }

    {
        using cell_index = unsigned int;
        int last_idx = 0;
        for (int _i = 0; _i < data_size; ++_i) {
            const auto &pos = positions[indices[_i]];
            const cell_index i = static_cast<const cell_index>(floor(pos.data[0] / .001));
            const cell_index j = static_cast<const cell_index>(floor(pos.data[1] / .001));
            const cell_index k = static_cast<const cell_index>(floor(pos.data[2] / .001));
            bitmask_t coords[3]{i, j, k};
            auto idx = static_cast<unsigned int>(hilbert_c2i(3, CHAR_BIT, coords));
            ASSERT_GE(idx, last_idx);
            last_idx = idx;
        }
    }

    {
        std::vector<std::size_t> inverseIndices(indices.size());
        for(std::size_t i = 0; i < indices.size(); ++i) {
            inverseIndices[indices[i]] = i;
        }
        readdy::util::collections::reorder_destructive(inverseIndices.begin(), inverseIndices.end(), positions.begin());
    }

    {
        using cell_index = unsigned int;
        int last_idx = 0;
        for (int _i = 0; _i < data_size; ++_i) {
            const auto &pos = positions[_i];
            const cell_index i = static_cast<const cell_index>(floor(pos.data[0] / .001));
            const cell_index j = static_cast<const cell_index>(floor(pos.data[1] / .001));
            const cell_index k = static_cast<const cell_index>(floor(pos.data[2] / .001));
            bitmask_t coords[3]{i, j, k};
            auto idx = static_cast<unsigned int>(hilbert_c2i(3, CHAR_BIT, coords));
            ASSERT_GE(idx, last_idx);
        }
    }

}

}