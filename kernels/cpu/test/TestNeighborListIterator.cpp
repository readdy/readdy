/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file TestNeighborListIterator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 15.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include "gtest/gtest.h"
#include <readdy/kernel/cpu/nl/CellLinkedList.h>

namespace {

TEST(TestNeighborListIterator, BoxIterator) {
    using namespace readdy;
    model::Context context;

    kernel::cpu::thread_pool pool (readdy_default_n_threads());
    kernel::cpu::data::DefaultDataContainer data (context, pool);
    kernel::cpu::nl::CompactCellLinkedList ccll(data, context, pool);

    context.particleTypes().add("A", 1.0);
    context.reactions().addFusion("fusion", "A", "A", "A", 1.0, 1.0);
    context.boxSize() = {5., 5., 5.};

    std::size_t n = 5;

    {
        std::vector<model::Particle> particles;
        for(auto i = 0; i < n; ++i) {
            particles.emplace_back(0, 0, 0, context.particleTypes().idOf("A"));
        }
        data.addParticles(particles);
    }

    std::vector<model::Particle::id_type> ids;
    for(std::size_t i=0; i < n; ++i){
        ids.push_back(data.getParticle(i).getId());
    }

    ccll.setUp(0, 1, {});

    {
        std::size_t nNonemptyCells {0};
        for(const auto &head : ccll.head()) {
            if((*head).load() != 0) {
                ++nNonemptyCells;
            }
        }
        ASSERT_EQ(nNonemptyCells, 1);
    }

    auto it = std::find_if(ccll.head().begin(), ccll.head().end(), [](const auto& head) {
        return (*head).load() != 0;
    });
    ASSERT_NE(it, ccll.head().end());

    auto cell = static_cast<std::size_t>(std::distance(ccll.head().begin(), it));
    for(auto boxIt = ccll.particlesBegin(cell); boxIt != ccll.particlesEnd(cell); ++boxIt) {
        auto idIt = std::find(ids.begin(), ids.end(), data.entry_at(*boxIt).id);
        ASSERT_NE(idIt, ids.end());
        ids.erase(idIt);
    }

    ASSERT_EQ(ids.size(), 0);
}

TEST(TestNeighborListIterator, BoxIteratorEmptyBox) {
    using namespace readdy;
    model::Context context;
    kernel::cpu::thread_pool pool (readdy_default_n_threads());
    kernel::cpu::data::DefaultDataContainer data (context, pool);
    kernel::cpu::nl::CompactCellLinkedList ccll(data, context, pool);

    context.particleTypes().add("A", 1.0);
    context.reactions().addFusion("fusion", "A", "A", "A", 1.0, 1.0);
    context.boxSize() = {5., 5., 5.};

    ccll.setUp(0, 1, {});
    ccll.update({});

    {
        std::size_t nNonemptyCells {0};
        for(const auto &head : ccll.head()) {
            if((*head).load() != 0) {
                ++nNonemptyCells;
            }
        }
        ASSERT_EQ(nNonemptyCells, 0);
    }

    // some random cell
    auto cell = 1_z;
    ASSERT_EQ(ccll.particlesBegin(cell), ccll.particlesEnd(cell));
}


}