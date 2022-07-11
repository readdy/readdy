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
 * @file TestCellLinkedList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12.09.17
 * @copyright BSD-3
 */

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <readdy/kernel/cpu/nl/ContiguousCellLinkedList.h>
#include <readdy/kernel/cpu/nl/CellLinkedList.h>

using namespace readdy;

using CompactCLL = readdy::kernel::cpu::nl::CompactCellLinkedList;
using ContiguousCLL = readdy::kernel::cpu::nl::ContiguousCellLinkedList;

TEMPLATE_TEST_CASE("Test cpu cell linked list", "[cpu]", CompactCLL, ContiguousCLL) {
    auto cllRadius = GENERATE(1, 2, 3);

    INFO(fmt::format("Testing with cll radius {}", cllRadius));

    model::Context context;
    context.particleTypes().add("Test", 1.);
    auto id = context.particleTypes().idOf("Test");
    scalar cutoff = 1;
    context.reactions().addFusion("Fusion", id, id, id, 1., cutoff);
    context.boxSize()[0] = 10;
    context.boxSize()[1] = 10;
    context.boxSize()[2] = 10;
    bool periodic = true;
    context.periodicBoundaryConditions()[0] = periodic;
    context.periodicBoundaryConditions()[1] = periodic;
    context.periodicBoundaryConditions()[2] = periodic;

    kernel::cpu::thread_pool pool (readdy_default_n_threads());

    kernel::cpu::data::DefaultDataContainer data (context, pool);

    auto nParticles = 1000;
    for(int i = 0; i < nParticles; ++i) {
        model::Particle particle(model::rnd::uniform_real<scalar>(-5, 5),
                                 model::rnd::uniform_real<scalar>(-5, 5),
                                 model::rnd::uniform_real<scalar>(-5, 5), id);
        data.addParticle(particle);
    }

    SECTION("Insertion") {

        std::vector<ParticleId> particleIds;
        std::for_each(data.begin(), data.end(), [&] (const auto &entry){ particleIds.push_back(entry.id); });

        std::unique_ptr<TestType> cll = std::make_unique<TestType>(data, context, pool);
        cll->setUp(context.calculateMaxCutoff(), static_cast<kernel::cpu::nl::CellLinkedList::cell_radius_type>(cllRadius));
        cll->update();

        {
            std::size_t pidx{0};
            for (const auto &e : data) {
                if (!e.deactivated) {
                    auto cell = cll->cellOfParticle(pidx);
                    auto foundIt = std::find(cll->particlesBegin(cell), cll->particlesEnd(cell), pidx);
                    // particle pidx should be found in cell
                    REQUIRE(foundIt != cll->particlesEnd(cell));
                }
                ++pidx;
            }
        }

        for (std::size_t cell = 0; cell < cll->nCells(); ++cell) {
            for(auto itParticle = cll->particlesBegin(cell); itParticle != cll->particlesEnd(cell); ++itParticle) {
                const auto &entry = cll->data().entry_at(*itParticle);

                auto foundIt = std::find(particleIds.begin(), particleIds.end(), entry.id);
                REQUIRE(foundIt != particleIds.end());
                particleIds.erase(foundIt);

                std::vector<std::size_t> neighbors;
                cll->forEachNeighbor(*itParticle, [&neighbors](auto neighbor) {
                    neighbors.push_back(neighbor);
                });

                std::size_t pidx = 0;
                for (const auto &e : data) {
                    if (pidx != *itParticle && bcs::dist(entry.pos, e.pos, context.boxSize(), context.periodicBoundaryConditions()) < cutoff) {
                        REQUIRE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                    }
                    ++pidx;
                }

            }
        }

        REQUIRE(particleIds.empty());
    }

    SECTION("Insertion and deactivation") {

        for(std::size_t i = 0; i < nParticles; ++i) {
            if(model::rnd::uniform_int(0, 1) == 0) {
                data.removeEntry(i);
            }
        }

        std::unique_ptr<TestType> cll = std::make_unique<TestType>(data, context, pool);
        cll->setUp(context.calculateMaxCutoff(), static_cast<kernel::cpu::nl::CellLinkedList::cell_radius_type>(cllRadius));
        cll->update();

        for (std::size_t cell = 0; cell < cll->nCells(); ++cell) {
            for (auto itParticle = cll->particlesBegin(cell); itParticle != cll->particlesEnd(cell); ++itParticle) {
                const auto &entry = cll->data().entry_at(*itParticle);
                REQUIRE_FALSE(entry.deactivated);
                std::vector<std::size_t> neighbors;

                cll->forEachNeighbor(*itParticle, [&](auto neighbor) {
                    const auto &neighborEntry = cll->data().entry_at(neighbor);
                    REQUIRE_FALSE(neighborEntry.deactivated);
                    neighbors.push_back(neighbor);
                });

                // check for every particle that is closer than cutoff*cutoff that it is in "neighbors" vector
                std::size_t pidx = 0;
                for (const auto &e : cll->data()) {
                    if (!e.deactivated) {
                        if (pidx != *itParticle && bcs::dist(entry.pos, e.pos, context.boxSize(), context.periodicBoundaryConditions()) < cutoff) {
                            REQUIRE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                        }
                    }
                    ++pidx;
                }

            }
        }
    }

    SECTION("Diffuse") {
        std::unique_ptr<TestType> cll = std::make_unique<TestType>(data, context, pool);
        cll->setUp(context.calculateMaxCutoff(), static_cast<kernel::cpu::nl::CellLinkedList::cell_radius_type>(cllRadius));
        cll->update();

        std::size_t n_steps = 3;
        for(std::size_t t = 0; t < n_steps; ++t) {
            for (std::size_t i = 0; i < nParticles; ++i) {
                cll->data().displace(i, 2. * model::rnd::normal3<readdy::scalar>(0, 1));
            }

            cll->update();


            for(std::size_t cell = 0; cell < cll->nCells(); ++cell) {
                for(auto itParticle = cll->particlesBegin(cell); itParticle != cll->particlesEnd(cell); ++itParticle) {
                    const auto &entry = cll->data().entry_at(*itParticle);
                    REQUIRE_FALSE(entry.deactivated);
                    std::vector<std::size_t> neighbors;
                    cll->forEachNeighbor(*itParticle, [&](auto neighborIdx) {
                        const auto &neighborEntry = cll->data().entry_at(neighborIdx);
                        REQUIRE_FALSE(neighborEntry.deactivated);
                        neighbors.push_back(neighborIdx);
                    });

                    // check for every particle that is closer than cutoff*cutoff that it is in "neighbors" vector
                    std::size_t pidx = 0;
                    for (const auto &e : cll->data()) {
                        if (!e.deactivated) {
                            if (pidx != *itParticle && bcs::dist(entry.pos, e.pos, context.boxSize(),
                                    context.periodicBoundaryConditions()) < cutoff) {
                                REQUIRE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                            }
                        }
                        ++pidx;
                    }

                }
            }
        }
    }
}
