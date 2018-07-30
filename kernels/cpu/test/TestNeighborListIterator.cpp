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
 * @file TestNeighborListIterator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 15.09.17
 * @copyright GPL-3
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
        ids.push_back(data.getParticle(i).id());
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