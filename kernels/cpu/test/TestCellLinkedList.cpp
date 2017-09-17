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
 * @file TestCellLinkedList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <gtest/gtest.h>

#include "readdy/kernel/cpu/nl/CLLNeighborList.h"

namespace {

class TestCLLImpl : public ::testing::TestWithParam<const char*> {};

TEST_P(TestCLLImpl, TestInsert) {
    using namespace readdy;

    model::KernelContext context;
    context.particle_types().add("Test", 1., 1.);
    auto id = context.particle_types().id_of("Test");
    scalar cutoff = 1;
    context.reactions().add(std::make_unique<model::reactions::Fusion>("Fusion", id, id, id, 1., cutoff));
    context.boxSize()[0] = 10;
    context.boxSize()[1] = 10;
    context.boxSize()[2] = 10;
    bool periodic = true;
    context.periodicBoundaryConditions()[0] = periodic;
    context.periodicBoundaryConditions()[1] = periodic;
    context.periodicBoundaryConditions()[2] = periodic;

    util::thread::Config config;
    config.setNThreads(1);
    
    kernel::cpu::nl::NeighborList *nl;
    
    kernel::cpu::nl::CompactCLLNeighborList ccll_nl (1, context, config);
    kernel::cpu::nl::DynamicCLLNeighborList dcll_nl (context, config);
    kernel::cpu::nl::ContiguousCLLNeighborList contiguous_cll_nl (1, context, config);
    
    if(GetParam() == "DynamicCLL") {
        nl = &dcll_nl;
    } else if(GetParam() == "CompactCLL") {
        nl = &ccll_nl;
    } else if(GetParam() == "ContiguousCLL") {
        nl = &contiguous_cll_nl;
    } else {
        throw std::invalid_argument("This should not happen");
    }

    for(int i= 0; i < 1000; ++i) {
        model::Particle particle(model::rnd::uniform_real<scalar>(-5, 5), 
                                 model::rnd::uniform_real<scalar>(-5, 5), 
                                 model::rnd::uniform_real<scalar>(-5, 5), id);
        nl->data()->addParticle(particle);
    }

    context.configure(false);

    nl->set_up({});
    nl->update({});
    
    const auto &d2 = context.distSquaredFun();

    {
        for (auto it = nl->begin(); it != nl->end(); ++it) {
            const auto &entry = nl->data()->entry_at(it->current_particle());
            std::vector<std::size_t> neighbors;
            neighbors.reserve(it->n_neighbors());
            for (unsigned long neighbor : *it) {
                const auto &neighborEntry = nl->data()->entry_at(neighbor);
                neighbors.push_back(neighbor);
                ASSERT_LE(d2(entry.pos, neighborEntry.pos), cutoff*cutoff);
            }

            // check for every particle that is closer than cutoff*cutoff that it is in "neighbors" vector
            std::size_t pidx = 0;
            for(const auto &e : *nl->data()) {
                if(pidx != it->current_particle() && d2(entry.pos, e.pos) < cutoff*cutoff) {
                    ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                }
                ++pidx;
            }

        }
    }
}

}

INSTANTIATE_TEST_CASE_P(TestCellLinkedList, TestCLLImpl,
                        ::testing::Values("DynamicCLL", "CompactCLL", "ContiguousCLL"));
