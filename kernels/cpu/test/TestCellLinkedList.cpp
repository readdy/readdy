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
#include <readdy/kernel/cpu/nl/ContiguousCellLinkedList.h>

#include "readdy/kernel/cpu/nl/CellLinkedList.h"

namespace {

class TestCLL : public ::testing::TestWithParam<::testing::tuple<int, const char*>> {

public:
    auto cllRadius() const {
        return static_cast<readdy::kernel::cpu::nl::CellLinkedList::cell_radius_type>(::testing::get<0>(GetParam()));
    }

    std::string cllName() const {
        return ::testing::get<1>(GetParam());
    }
};

TEST_P(TestCLL, Insert) {
    using namespace readdy;

    model::Context context;
    context.particle_types().add("Test", 1.);
    auto id = context.particle_types().idOf("Test");
    scalar cutoff = 1;
    context.reactions().addFusion("Fusion", id, id, id, 1., cutoff);
    context.boxSize()[0] = 10;
    context.boxSize()[1] = 10;
    context.boxSize()[2] = 10;
    bool periodic = true;
    context.periodicBoundaryConditions()[0] = periodic;
    context.periodicBoundaryConditions()[1] = periodic;
    context.periodicBoundaryConditions()[2] = periodic;

    util::thread::Config config;

    kernel::cpu::data::DefaultDataContainer data (context, config);

    for(int i = 0; i < 1000; ++i) {
        model::Particle particle(model::rnd::uniform_real<scalar>(-5, 5), 
                                 model::rnd::uniform_real<scalar>(-5, 5), 
                                 model::rnd::uniform_real<scalar>(-5, 5), id);
        data.addParticle(particle);
    }

    std::vector<model::Particle::id_type> particleIds;
    std::for_each(data.begin(), data.end(), [&] (const auto &entry){ particleIds.push_back(entry.id); });

    context.configure();
    const auto &d2 = context.distSquaredFun();

    if(cllName() == "CompactCLL") {
        kernel::cpu::nl::CompactCellLinkedList ccll(data, context, config);
        ccll.setUp(0, cllRadius(), {});
        ccll.update({});

        for (std::size_t cell = 0; cell < ccll.nCells(); ++cell) {
            for (auto it = ccll.cellNeighborsBegin(cell); it != ccll.cellNeighborsEnd(cell); ++it) {
                const auto &entry = ccll.data().entry_at(it.currentParticle());
                auto foundIt = std::find(particleIds.begin(), particleIds.end(), entry.id);
                ASSERT_NE(foundIt, particleIds.end());
                particleIds.erase(foundIt);
                std::vector<std::size_t> neighbors;
                for (auto itNeigh = it.neighborsBegin(); itNeigh != it.neighborsEnd(); ++itNeigh) {
                    const auto &neighborEntry = ccll.data().entry_at(*itNeigh);
                    neighbors.push_back(*itNeigh);
                }

                // check for every particle that is closer than cutoff*cutoff that it is in "neighbors" vector
                std::size_t pidx = 0;
                for (const auto &e : ccll.data()) {
                    if (pidx != it.currentParticle() && d2(entry.pos, e.pos) < cutoff * cutoff) {
                        ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                    }
                    ++pidx;
                }

            }
        }
    } else if(cllName() == "ContiguousCLL") {
        kernel::cpu::nl::ContiguousCellLinkedList ccll(data, context, config);
        ccll.setUp(0, cllRadius(), {});
        ccll.update({});

        {
            std::size_t pidx{0};
            for (const auto &e : data) {
                if (!e.deactivated) {
                    auto cell = ccll.cellOfParticle(pidx);
                    auto foundIt = std::find(ccll.particlesBegin(cell), ccll.particlesEnd(cell), pidx);
                    ASSERT_NE(foundIt, ccll.particlesEnd(cell)) << "Did not find particle " << pidx
                                                                << " in cell " << cell;
                }
                ++pidx;
            }
        }

        for (std::size_t cell = 0; cell < ccll.nCells(); ++cell) {
            for(auto itParticle = ccll.particlesBegin(cell); itParticle != ccll.particlesEnd(cell); ++itParticle) {
                const auto &entry = ccll.data().entry_at(*itParticle);

                auto foundIt = std::find(particleIds.begin(), particleIds.end(), entry.id);
                ASSERT_NE(foundIt, particleIds.end());
                particleIds.erase(foundIt);

                std::vector<std::size_t> neighbors;
                ccll.forEachNeighbor(*itParticle, [&neighbors](auto neighbor) {
                    neighbors.push_back(neighbor);
                });

                std::size_t pidx = 0;
                for (const auto &e : data) {
                    if (pidx != *itParticle && d2(entry.pos, e.pos) < cutoff * cutoff) {
                        ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                    }
                    ++pidx;
                }

            }
        }

    }
    ASSERT_TRUE(particleIds.empty());
}

TEST_P(TestCLL, InsertAndDeactivate) {
    using namespace readdy;

    model::Context context;
    context.particle_types().add("Test", 1.);
    auto id = context.particle_types().idOf("Test");
    scalar cutoff = 1;
    context.reactions().addFusion("Fusion", id, id, id, 1., cutoff);
    context.boxSize()[0] = 10;
    context.boxSize()[1] = 10;
    context.boxSize()[2] = 10;
    bool periodic = true;
    context.periodicBoundaryConditions()[0] = periodic;
    context.periodicBoundaryConditions()[1] = periodic;
    context.periodicBoundaryConditions()[2] = periodic;

    util::thread::Config config;

    kernel::cpu::data::DefaultDataContainer data (context, config);

    auto n_particles = 1000;
    for(int i = 0; i < n_particles; ++i) {
        model::Particle particle(model::rnd::uniform_real<scalar>(-5, 5),
                                 model::rnd::uniform_real<scalar>(-5, 5),
                                 model::rnd::uniform_real<scalar>(-5, 5), id);
        data.addParticle(particle);
    }


    for(std::size_t i = 0; i < n_particles; ++i) {
        if(model::rnd::uniform_int(0, 1) == 0) {
            data.removeEntry(i);
        }
    }

    context.configure();
    const auto &d2 = context.distSquaredFun();

    if(cllName() == "CompactCLL") {
        kernel::cpu::nl::CompactCellLinkedList nl(data, context, config);
        nl.setUp(0, cllRadius(), {});
        nl.update({});

        for (std::size_t cell = 0; cell < nl.nCells(); ++cell) {
            for (auto it = nl.cellNeighborsBegin(cell); it != nl.cellNeighborsEnd(cell); ++it) {
                const auto &entry = nl.data().entry_at(it.currentParticle());
                ASSERT_FALSE(entry.deactivated);
                std::vector<std::size_t> neighbors;
                for (auto itNeigh = it.neighborsBegin(); itNeigh != it.neighborsEnd(); ++itNeigh) {
                    const auto &neighborEntry = nl.data().entry_at(*itNeigh);
                    ASSERT_FALSE(neighborEntry.deactivated);
                    neighbors.push_back(*itNeigh);
                }

                // check for every particle that is closer than cutoff*cutoff that it is in "neighbors" vector
                std::size_t pidx = 0;
                for (const auto &e : nl.data()) {
                    if (!e.deactivated) {
                        if (pidx != it.currentParticle() && d2(entry.pos, e.pos) < cutoff * cutoff) {
                            ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                        }
                    }
                    ++pidx;
                }

            }
        }
    } else if(cllName() == "ContiguousCLL") {
        kernel::cpu::nl::ContiguousCellLinkedList nl(data, context, config);
        nl.setUp(0, cllRadius(), {});
        nl.update({});

        for (std::size_t cell = 0; cell < nl.nCells(); ++cell) {
            for (auto itParticle = nl.particlesBegin(cell); itParticle != nl.particlesEnd(cell); ++itParticle) {
                const auto &entry = nl.data().entry_at(*itParticle);
                ASSERT_FALSE(entry.deactivated);
                std::vector<std::size_t> neighbors;

                nl.forEachNeighbor(*itParticle, [&](auto neighbor) {
                    const auto &neighborEntry = nl.data().entry_at(neighbor);
                    ASSERT_FALSE(neighborEntry.deactivated);
                    neighbors.push_back(neighbor);
                });

                // check for every particle that is closer than cutoff*cutoff that it is in "neighbors" vector
                std::size_t pidx = 0;
                for (const auto &e : nl.data()) {
                    if (!e.deactivated) {
                        if (pidx != *itParticle && d2(entry.pos, e.pos) < cutoff * cutoff) {
                            ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                        }
                    }
                    ++pidx;
                }

            }
        }
    }
}

TEST_P(TestCLL, Diffuse) {
    using namespace readdy;

    model::Context context;
    context.particle_types().add("Test", 1.);
    auto id = context.particle_types().idOf("Test");
    scalar cutoff = 1;
    context.reactions().addFusion("Fusion", id, id, id, 1., cutoff);
    context.boxSize()[0] = 10;
    context.boxSize()[1] = 10;
    context.boxSize()[2] = 10;
    bool periodic = true;
    context.periodicBoundaryConditions()[0] = periodic;
    context.periodicBoundaryConditions()[1] = periodic;
    context.periodicBoundaryConditions()[2] = periodic;

    util::thread::Config config;

    kernel::cpu::data::DefaultDataContainer data (context, config);

    auto n_particles = 1000;
    for(int i = 0; i < n_particles; ++i) {
        model::Particle particle(model::rnd::uniform_real<scalar>(-5, 5),
                                 model::rnd::uniform_real<scalar>(-5, 5),
                                 model::rnd::uniform_real<scalar>(-5, 5), id);
        data.addParticle(particle);
    }

    context.configure();

    if (cllName() == "CompactCLL") {
        kernel::cpu::nl::CompactCellLinkedList nl(data, context, config);
        nl.setUp(0, cllRadius(), {});
        nl.update({});

        std::size_t n_steps = 3;
        for(std::size_t t = 0; t < n_steps; ++t) {
            for (std::size_t i = 0; i < n_particles; ++i) {
                nl.data().displace(i, 2. * model::rnd::normal3<readdy::scalar>(0, 1));
            }

            nl.update({});

            const auto &d2 = context.distSquaredFun();

            for(std::size_t cell = 0; cell < nl.nCells(); ++cell) {
                for (auto it = nl.cellNeighborsBegin(cell); it != nl.cellNeighborsEnd(cell); ++it) {
                    const auto &entry = nl.data().entry_at(it.currentParticle());
                    ASSERT_FALSE(entry.deactivated);
                    std::vector<std::size_t> neighbors;
                    for (auto itNeigh = it.neighborsBegin(); itNeigh != it.neighborsEnd(); ++itNeigh) {
                        const auto &neighborEntry = nl.data().entry_at(*itNeigh);
                        ASSERT_FALSE(neighborEntry.deactivated);
                        neighbors.push_back(*itNeigh);
                    }

                    // check for every particle that is closer than cutoff*cutoff that it is in "neighbors" vector
                    std::size_t pidx = 0;
                    for (const auto &e : nl.data()) {
                        if (!e.deactivated) {
                            if (pidx != it.currentParticle() && d2(entry.pos, e.pos) < cutoff * cutoff) {
                                ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                            }
                        }
                        ++pidx;
                    }

                }
            }
        }
    } else if(cllName() == "ContiguousCLL") {
        kernel::cpu::nl::ContiguousCellLinkedList nl(data, context, config);
        nl.setUp(0, cllRadius(), {});
        nl.update({});

        std::size_t n_steps = 3;
        for(std::size_t t = 0; t < n_steps; ++t) {
            for (std::size_t i = 0; i < n_particles; ++i) {
                nl.data().displace(i, 2. * model::rnd::normal3<readdy::scalar>(0, 1));
            }

            nl.update({});

            const auto &d2 = context.distSquaredFun();

            for(std::size_t cell = 0; cell < nl.nCells(); ++cell) {
                for(auto itParticle = nl.particlesBegin(cell); itParticle != nl.particlesEnd(cell); ++itParticle) {
                    const auto &entry = nl.data().entry_at(*itParticle);
                    ASSERT_FALSE(entry.deactivated);
                    std::vector<std::size_t> neighbors;
                    nl.forEachNeighbor(*itParticle, [&](auto neighborIdx) {
                        const auto &neighborEntry = nl.data().entry_at(neighborIdx);
                        ASSERT_FALSE(neighborEntry.deactivated);
                        neighbors.push_back(neighborIdx);
                    });

                    // check for every particle that is closer than cutoff*cutoff that it is in "neighbors" vector
                    std::size_t pidx = 0;
                    for (const auto &e : nl.data()) {
                        if (!e.deactivated) {
                            if (pidx != *itParticle && d2(entry.pos, e.pos) < cutoff * cutoff) {
                                ASSERT_NE(std::find(neighbors.begin(), neighbors.end(), pidx), neighbors.end());
                            }
                        }
                        ++pidx;
                    }

                }
            }
        }
    }

}

}

INSTANTIATE_TEST_CASE_P(TestCellLinkedList, TestCLL,
                        ::testing::Combine(::testing::Values(1,2,3), ::testing::Values("CompactCLL", "ContiguousCLL")));
