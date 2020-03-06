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


#include <array>

#include <catch2/catch.hpp>

#include <readdy/testing/Utils.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/common/numeric.h>

using namespace readdy;
using namespace readdytesting::kernel;

/**
 * << detailed description >>
 *
 * @file TestTopologyGraphs.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.03.17
 * @copyright BSD-3
 */

template<typename Tup>
bool containsTuple(const std::vector<Tup> &v, const Tup &t) {
    auto tReverse = readdy::util::reverse(t);
    auto it = std::find_if(std::begin(v), std::end(v), [&](const Tup &other) {
        return t == other || tReverse == other;
    });
    return it != std::end(v);
}

TEST_CASE("Test topology graphs") {

    SECTION("Test the hasher") {
        SECTION("Quadruples") {
            readdy::util::particle_type_quadruple_hasher quadHasher;
            std::array<readdy::ParticleTypeId, 4> range{1, 2, 3, 4};
            do {
                REQUIRE(util::sortTypeQuadruple(range[0], range[1], range[2], range[3]) == std::make_tuple(1, 2, 3, 4));
            } while (std::next_permutation(range.begin(), range.end()));
            REQUIRE(quadHasher(std::make_tuple(1, 2, 3, 4)) == quadHasher(std::make_tuple(4, 3, 2, 1)));
            REQUIRE(quadHasher(std::make_tuple(1, 3, 2, 4)) == quadHasher(std::make_tuple(4, 2, 3, 1)));
        }
        SECTION("Triples") {
            readdy::util::particle_type_triple_hasher tripleHasher;
            std::array<readdy::ParticleTypeId, 3> range{1, 2, 3};
            do {
                REQUIRE(readdy::util::sortTypeTriple(range[0], range[1], range[2]) == std::make_tuple(1, 2, 3));
            } while (std::next_permutation(range.begin(), range.end()));
            REQUIRE(tripleHasher(std::make_tuple(1, 2, 3)) == tripleHasher(std::make_tuple(3, 2, 1)));
            REQUIRE(tripleHasher(std::make_tuple(2, 1, 3)) == tripleHasher(std::make_tuple(3, 1, 2)));
        }
        SECTION("Tuples") {
            readdy::util::particle_type_pair_hasher pairHasher;
            REQUIRE(pairHasher(std::make_tuple(1, 2)) == pairHasher(std::make_tuple(2, 1)));
            util::particle_type_pair_unordered_map<int> map;
            int a = 1, b = 2;
            map[std::make_tuple(a, b)] = 5;
            REQUIRE(map[std::make_tuple(1, 2)] == 5);
            REQUIRE(map[std::make_tuple(2, 1)] == 5);
        }
    }
    SECTION("Test the graph") {
        readdy::model::top::Graph graph;
        auto it1 = graph.addVertex({{0, ""}});
        auto it2 = graph.addVertex({{1, ""}});
        SECTION("Access with indices") {
            graph.addEdge(it1, it2);
            REQUIRE(graph.vertices().size() == 2);
            REQUIRE(graph.vertices().at(it1).data().particleIndex == 0);
            REQUIRE(graph.vertices().at(it2).data().particleIndex == 1);
            REQUIRE(graph.vertices().at(it1).neighbors().size() == 1);
            REQUIRE(graph.vertices().at(it2).neighbors().size() == 1);
            REQUIRE(graph.vertices().at(graph.vertices().at(it1).neighbors()[0])->particleIndex == 1);
            REQUIRE(graph.vertices().at(graph.vertices().at(it1).neighbors()[0])->particleIndex == 0);
            graph.removeVertex(it1);
            REQUIRE(graph.vertices().size() == 1);
            REQUIRE(graph.begin()->data().particleIndex == 1);
            REQUIRE(graph.begin()->neighbors().empty());
        }
        SECTION("Connected components") {
            graph.addVertex({{2, ""}});

            auto vertex_ref_0 = graph.vertices().begin();
            auto vertex_ref_1 = ++graph.vertices().begin();
            auto vertex_ref_2 = ++(++graph.vertices().begin());

            auto it = graph.vertices().begin();
            auto it_adv = ++graph.vertices().begin();
            graph.addEdge(it, it_adv);

            auto subGraphs = graph.connectedComponents();
            REQUIRE(subGraphs.size() == 2);
            {
                REQUIRE(subGraphs[0].vertices().size() == 2);
                REQUIRE(subGraphs[0].vertices().begin() == vertex_ref_0);
                REQUIRE(++subGraphs[0].vertices().begin() == vertex_ref_1);
            }
            {
                REQUIRE(subGraphs[1].vertices().size() == 1);
                REQUIRE(subGraphs[1].vertices().begin() == vertex_ref_2);
            }
        }
        SECTION("Find N tuples") {
            graph.addVertex({{2, ""}});
            graph.addVertex({{3, ""}});

            auto it = graph.begin();
            graph.addEdge(it, ++it);
            graph.addEdge(it, ++it);
            graph.addEdge(it, ++it);
            graph.addEdge(it, graph.begin());

            auto n_tuples = graph.findNTuples();
            const auto &pairs = std::get<0>(n_tuples);
            const auto &triples = std::get<1>(n_tuples);
            const auto &quadruples = std::get<2>(n_tuples);

            REQUIRE(pairs.size() == 4);
            REQUIRE(triples.size() == 4);
            REQUIRE(quadruples.size() == 4);

            auto a = graph.begin().persistent_index();
            auto b = std::next(graph.begin()).persistent_index();
            auto c = std::next(graph.begin(), 2).persistent_index();
            auto d = std::next(graph.begin(), 3).persistent_index();

            {
                // tuples
                REQUIRE(containsTuple(pairs, std::make_tuple(a, b)));
                REQUIRE(containsTuple(pairs, std::make_tuple(b, c)));
                REQUIRE(containsTuple(pairs, std::make_tuple(c, d)));
                REQUIRE(containsTuple(pairs, std::make_tuple(d, a)));
            }

            {
                // triples
                REQUIRE(containsTuple(triples, std::make_tuple(a, b, c)));
                REQUIRE(containsTuple(triples, std::make_tuple(b, c, d)));
                REQUIRE(containsTuple(triples, std::make_tuple(c, d, a)));
                REQUIRE(containsTuple(triples, std::make_tuple(d, a, b)));
            }

            {
                // quadruple
                REQUIRE(containsTuple(quadruples, std::make_tuple(d, a, b, c)));
                REQUIRE(containsTuple(quadruples, std::make_tuple(a, b, c, d)));
                REQUIRE(containsTuple(quadruples, std::make_tuple(b, c, d, a)));
                REQUIRE(containsTuple(quadruples, std::make_tuple(c, d, a, b)));
            }
        }
        SECTION("Find N tuples in triangle") {
            graph.addVertex({{2}});

            graph.addEdge(graph.begin(), std::next(graph.begin()));
            graph.addEdge(std::next(graph.begin()), std::next(graph.begin(), 2));
            graph.addEdge(std::next(graph.begin(), 2), graph.begin());

            auto n_tuples = graph.findNTuples();
            const auto &tuples = std::get<0>(n_tuples);
            const auto &triples = std::get<1>(n_tuples);
            const auto &quadruples = std::get<2>(n_tuples);

            REQUIRE(tuples.size() == 3);
            REQUIRE(triples.size() == 3);
            REQUIRE(quadruples.empty());

            auto a = graph.begin().persistent_index();
            auto b = std::next(graph.begin()).persistent_index();
            auto c = std::next(std::next(graph.begin())).persistent_index();

            // tuples
            REQUIRE(containsTuple(tuples, std::make_tuple(a, b)));
            REQUIRE(containsTuple(tuples, std::make_tuple(b, c)));
            REQUIRE(containsTuple(tuples, std::make_tuple(c, a)));

            // triples
            REQUIRE(containsTuple(triples, std::make_tuple(a, b, c)));
            REQUIRE(containsTuple(triples, std::make_tuple(b, c, a)));
            REQUIRE(containsTuple(triples, std::make_tuple(c, a, b)));
        }
        SECTION("Are connected with N or less edges") {
            for (unsigned i = 2; i < 7; i++) {
                graph.addVertex({{i}});
            }
            // setup linear graph
            for (unsigned i = 0; i < 6; i++) {
                graph.addEdge(std::next(graph.begin(), i), std::next(graph.begin(), i + 1));
            }
            for (unsigned i = 1; i < 7; i++) {
                auto d = graph.graphDistance(graph.begin(), std::next(graph.begin(), i));
                REQUIRE(d != -1);
                REQUIRE(d == i);
            }
            // now connect both ends of the graph
            graph.addEdge(std::next(graph.begin(), 6), graph.begin());
            for (unsigned i = 0; i < 5; i++) {
                REQUIRE(graph.graphDistance(std::next(graph.begin(), i), std::next(graph.begin(), i+2)) == 2);
                REQUIRE(graph.graphDistance(std::next(graph.begin(), i), std::next(graph.begin(), (i-2)%7)) == 2);
                /*REQUIRE(graph.areConnectedWithNOrLessEdges(2, *(std::next(graph.firstVertex(), i)),
                                                           *std::next(graph.firstVertex(), i + 2)));
                REQUIRE(graph.areConnectedWithNOrLessEdges(2, *std::next(graph.firstVertex(), i),
                                                           *std::next(graph.firstVertex(), (i - 2) % 7)));
                REQUIRE(!graph.areConnectedWithNOrLessEdges(1, *(std::next(graph.firstVertex(), i)),
                                                            *std::next(graph.firstVertex(), i + 2)));
                REQUIRE(!graph.areConnectedWithNOrLessEdges(1, *std::next(graph.firstVertex(), i),
                                                            *std::next(graph.firstVertex(), (i - 2) % 7)));*/

            }
        }
    }

    SECTION("Append particle") {
        using namespace readdy;
        model::Context context;
        context.topologyRegistry().potentialConfiguration().pairPotentials[std::make_tuple(0, 0)].emplace_back();
        context.topologyRegistry().potentialConfiguration().pairPotentials[std::make_tuple(0, 1)].emplace_back();
        model::top::Graph graph;
        graph.addVertex({{10}});
        graph.addVertex({{1}});
        graph.addVertex({{200}});
        model::top::GraphTopology gt{0, graph, context, nullptr};
        {
            auto it = gt.graph().vertices().begin();
            auto it2 = ++gt.graph().vertices().begin();
            gt.addEdge(it.persistent_index(), it2.persistent_index());
            gt.addEdge((++it).persistent_index(), (++it2).persistent_index());
        }
        gt.configure();
        gt.appendParticle(13, 1);
        gt.configure();

        auto pred = [](const model::top::Vertex &v) -> bool {
            return v->particleIndex == 3;
        };
        auto it = std::find_if(gt.graph().begin(), gt.graph().end(), pred);
        REQUIRE(it != gt.graph().vertices().end());
        auto v2 = std::next(gt.graph().vertices().begin());
        REQUIRE(gt.containsEdge(it.persistent_index(), v2.persistent_index()));
    }

}
