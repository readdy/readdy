//
// Created by mho on 10/28/19.
//

#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include <catch2/catch.hpp>

#include <graphs/graphs.h>

auto debug = [](const std::string &str) {
    std::cerr << "DEBUG: " << str << std::endl;
};

template<typename T, size_t... I>
auto reverse_impl(T t, std::index_sequence<I...>) {
    return std::make_tuple(std::get<sizeof...(I) - 1 - I>(std::forward<T>(t))...);
}

template<typename T>
auto reverse(T t) {
    return reverse_impl(std::forward<T>(t), std::make_index_sequence<std::tuple_size<T>::value>());
}

template<typename T, size_t... I>
auto to_persistent_ix_impl(T t, std::index_sequence<I...>) {
    return std::make_tuple(graphs::PersistentIndex{static_cast<std::size_t>(std::get<I>(std::forward<T>(t)))}...);
}

template<typename T>
auto to_persistent_ix(T t) {
    return to_persistent_ix_impl(std::forward<T>(t), std::make_index_sequence<std::tuple_size<T>::value>());
}

template<typename T, size_t... I>
auto to_active_ix_impl(T t, std::index_sequence<I...>) {
    return std::make_tuple(std::get<I>(std::forward<T>(t)).value...);
}

template<typename T>
auto to_active_ix(T t) {
    return to_active_ix_impl(std::forward<T>(t), std::make_index_sequence<std::tuple_size<T>::value>());
}

std::size_t n_choose_k(std::size_t n, std::size_t k) {
    if(k == 0 || k == n) return 1;
    return n_choose_k(n - 1, k - 1) + n_choose_k(n - 1, k);
}

template<typename T1, typename T2>
std::size_t nTupleOccurrences(const std::vector<T1> &v, const T2 &t) {
    std::size_t n = 0;
    auto tReverse = reverse<T2>(t);
    for(auto it = std::begin(v); it != std::end(v); ++it) {
        if(*it == t || *it == tReverse) {
            ++n;
        }
    }
    return n;
}

template<typename T1, typename T2>
bool containsTupleXOR(const std::vector<T1> &v, const T2 &t) {
    return nTupleOccurrences(v, t) == 1;
}

graphs::DefaultGraph fullyConnectedGraph(std::size_t size) {
    graphs::DefaultGraph graph;
    for(std::size_t i = 0; i < size; ++i) {
        graph.addVertex(i);
    }
    for(std::size_t i = 0; i < graph.vertices().size(); ++i) {
        for(std::size_t j = i+1; j < graph.vertices().size(); ++j) {
            graph.addEdge(i, j);
        }
    }
    return graph;
}

template<std::size_t K>
auto quadruplesFullyConnectedPrimitive() {
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> quads;

    for(std::size_t i = 0; i < K; ++i){
        for(std::size_t j = 0; j < K; ++j){
            for(std::size_t k = 0; k < K; ++k){
                for(std::size_t l = 0; l < K; ++l){
                    std::array<std::size_t, 4> arr {{i, j, k, l}};
                    std::sort(std::begin(arr), std::end(arr));
                    auto pos = std::adjacent_find(std::begin(arr), std::end(arr));
                    if( pos == std::end(arr) ) {
                        // only true if no duplicates
                        auto tup = std::make_tuple(i, j, k, l);
                        if (!containsTupleXOR(quads, tup)) {
                            quads.push_back(tup);
                        }
                    }
                }
            }
        }
    }

    return quads;
}

template<std::size_t K>
void testFullyConnected() {
    GIVEN("A fully connected graph of size " + std::to_string(K)) {
        auto graph = fullyConnectedGraph(K);
        THEN("The number of vertices should be " + std::to_string(K)) {
            REQUIRE(graph.vertices().size() == K);
        }
        THEN("The number of neighbors per vertex should be " + std::to_string(K - 1)) {
            for(auto it = graph.vertices().begin(); it != graph.vertices().end(); ++it) {
                REQUIRE(it->neighbors().size() == K - 1);
            }
        }
        THEN("The number of edges should be 0.5n(n-1)") {
            REQUIRE(graph.nEdges() == static_cast<std::size_t>(K * (K - 1) / 2));
        }
        const auto &[pairs, triples, quadruples] = graph.findNTuples();
        THEN("The number of unique pairs should be 0.5n(n-1)") {
            REQUIRE(pairs.size() == static_cast<std::size_t>(K * (K - 1) / 2));
            for(std::size_t i = 0; i < graph.vertices().size(); ++i) {
                for (std::size_t j = i+1; j < graph.vertices().size(); ++j) {
                    REQUIRE(containsTupleXOR(pairs, std::make_tuple(graphs::DefaultGraph::PersistentVertexIndex{i},
                                                                    graphs::DefaultGraph::PersistentVertexIndex{j})));
                }
            }
        }
        THEN("The number of unique triples should be 3*(n choose 3)") {
            // Since the graph is fully connected (n choose 3) gives the number of different paths of length 3
            // in the graph. The number is in terms of sets, i.e., a path is defined via its vertices independent of
            // their order. For a set of vertices {v1, v2, v3} there are three different ways of enumeration
            //  - v1 v2 v3
            //  - v1 v3 v2
            //  - v2 v1 v3
            // which give distinct potential configurations.
            REQUIRE(triples.size() == 3*n_choose_k(K, 3));
            for(const auto &triple : triples) {
                CAPTURE(triple);
                REQUIRE(nTupleOccurrences(triples, triple) == 1);
            }
        }

        THEN("The number of unique quadruples should be 12*(n choose 4)") {
            // See number of unique triplets but this time there are 12 ways of enumeration:
            // (n choose 4) different paths of length 4, 4! * (n choose 4) different order dependent paths,
            // 0.5 * 4! * (n choose 4) = 12 * (n choose 4) paths where (v1 v2 v3 v4) == (v4 v3 v2 v1).
            auto quadsPrimitive = quadruplesFullyConnectedPrimitive<K>();
            REQUIRE(quadsPrimitive.size() == 12*n_choose_k(K, 4));
            REQUIRE(quadruples.size() == 12*n_choose_k(K, 4));

            for(const auto &quad : quadruples) {
                const auto& [v1, v2, v3, v4] = quad;
                CAPTURE(quad);
                REQUIRE(nTupleOccurrences(quadruples, quad) == 1);
                REQUIRE(nTupleOccurrences(quadsPrimitive, to_active_ix(std::make_tuple(v1, v2, v3, v4))) == 1);
            }
        }
    }
}

SCENARIO("Testing graphs basic functionality", "[graphs]") {

    GIVEN("Some tuples (A, B, C, D), (D, C, B, A), (A, C, B, D)") {
        auto t1 = std::make_tuple(0, 1, 2, 3);
        auto t2 = std::make_tuple(3, 2, 1, 0);
        auto t3 = std::make_tuple(0, 2, 1, 3);
        auto vec = std::vector{t1, t2, t3};
        THEN("the number of occurrences is 2, 2, 1") {
            REQUIRE(nTupleOccurrences(vec, t1) == 2);
            REQUIRE(nTupleOccurrences(vec, t2) == 2);
            REQUIRE(nTupleOccurrences(vec, t3) == 1);
        }
        THEN("containsTupleXOR is false, false, true") {
            REQUIRE_FALSE(containsTupleXOR(vec, t1));
            REQUIRE_FALSE(containsTupleXOR(vec, t2));
            REQUIRE(containsTupleXOR(vec, t3));
        }
    }

    GIVEN("A graph with source and target and two paths of different length between them") {
        graphs::Graph<graphs::IndexPersistentVector, graphs::Vertex<std::size_t>> graph;
        auto source = graph.addVertex(0);
        auto v1 = graph.addVertex(1);
        auto v2 = graph.addVertex(2);
        auto v3 = graph.addVertex(3);
        auto target = graph.addVertex(4);

        graph.addEdge(source, v1);
        graph.addEdge(v1, v2);
        graph.addEdge(v2, target);

        graph.addEdge(source, v3);
        graph.addEdge(v3, target);

        REQUIRE(graph.graphDistance(source, target) == 2);
    }

    GIVEN("A graph with two vertices") {
        graphs::Graph<graphs::IndexPersistentVector, graphs::Vertex<std::size_t>> graph;
        graph.addVertex(0);
        graph.addVertex(1);
        WHEN("connecting the two vertices") {
            graph.addEdge(0, 1);
            THEN("this should be reflected in the neighbors structure accessed by particle indices") {
                REQUIRE(graph.isConnected());
                REQUIRE(graph.vertices().size() == 2);
                REQUIRE(graph.vertices().begin()->neighbors().size() == 1);
                REQUIRE((++graph.vertices().begin())->neighbors().size() == 1);
                REQUIRE(graph.vertices().begin()->neighbors()[0].value == 1);
                REQUIRE((++graph.vertices().begin())->neighbors()[0].value == 0);
            }
            WHEN("removing the first particle") {
                graph.removeVertex(0);
                THEN("the size of the graph is 1") {
                    REQUIRE(graph.vertices().size() == 1);
                    REQUIRE(graph.vertices().begin()->neighbors().empty());
                    REQUIRE(graph.vertices().begin()->data() == 1);
                    REQUIRE(graph.nEdges() == 0);
                }
            }
        }
        WHEN("Adding a third vertex, connecting (0 -- 1), (2)") {
            graph.addVertex(2);

            graph.addEdge(0, 1);

            REQUIRE_FALSE(graph.isConnected());

            auto subGraphs = graph.connectedComponents();
            THEN("There should be two connected components") {
                REQUIRE(subGraphs.size() == 2);
            }
            THEN("The two connected components should be (0 -- 1) and (2)") {
                {
                    REQUIRE(subGraphs[0].vertices().size() == 2);
                    REQUIRE(subGraphs[0].vertices().begin()->data() == 0);
                    REQUIRE((++subGraphs[0].vertices().begin())->data() == 1);
                }
                {
                    REQUIRE(subGraphs[1].vertices().size() == 1);
                    REQUIRE(subGraphs[1].vertices().begin()->data() == 2);
                }
            }
        }
        WHEN("Adding two vertices and connecting (0 -- 1 -- 2 -- 3 -- 0)") {
            graph.addVertex(2);
            graph.addVertex(3);

            graph.addEdge(0, 1);
            graph.addEdge(1, 2);
            graph.addEdge(2, 3);
            graph.addEdge(3, 0);

            REQUIRE(graph.isConnected());

            auto n_tuples = graph.findNTuples();
            const auto& tuples = std::get<0>(n_tuples);
            const auto& triples = std::get<1>(n_tuples);
            const auto& quadruples = std::get<2>(n_tuples);

            THEN("expect 4 unique tuples") {
                REQUIRE(tuples.size() == 4);

                REQUIRE(containsTupleXOR(tuples, to_persistent_ix(std::make_tuple(0, 1))));
                REQUIRE(containsTupleXOR(tuples, to_persistent_ix(std::make_tuple(1, 2))));
                REQUIRE(containsTupleXOR(tuples, to_persistent_ix(std::make_tuple(2, 3))));
                REQUIRE(containsTupleXOR(tuples, to_persistent_ix(std::make_tuple(3, 0))));
            }
            THEN("expect 4 unique triples") {
                REQUIRE(triples.size() == 4);

                REQUIRE(containsTupleXOR(triples, to_persistent_ix(std::make_tuple(0, 1, 2))));
                REQUIRE(containsTupleXOR(triples, to_persistent_ix(std::make_tuple(1, 2, 3))));
                REQUIRE(containsTupleXOR(triples, to_persistent_ix(std::make_tuple(2, 3, 0))));
                REQUIRE(containsTupleXOR(triples, to_persistent_ix(std::make_tuple(3, 0, 1))));
            }
            THEN("expect 4 unique quadruples") {
                REQUIRE(quadruples.size() == 4);

                REQUIRE(containsTupleXOR(quadruples, to_persistent_ix(std::make_tuple(3, 0, 1, 2))));
                REQUIRE(containsTupleXOR(quadruples, to_persistent_ix(std::make_tuple(0, 1, 2, 3))));
                REQUIRE(containsTupleXOR(quadruples, to_persistent_ix(std::make_tuple(1, 2, 3, 0))));
                REQUIRE(containsTupleXOR(quadruples, to_persistent_ix(std::make_tuple(2, 3, 0, 1))));
            }
        }

        WHEN("Adding one vertex and connecting (0 -- 1 -- 2 -- 0)") {
            graph.addVertex(2);

            auto a = 0;
            auto b = 1;
            auto c = 2;

            graph.addEdge(a, b);
            graph.addEdge(b, c);
            graph.addEdge(c, a);

            auto n_tuples = graph.findNTuples();
            const auto& pairs = std::get<0>(n_tuples);
            const auto& triples = std::get<1>(n_tuples);
            const auto& quadruples = std::get<2>(n_tuples);

            REQUIRE(graph.isConnected());

            THEN("Expect 3 unique pairs") {
                REQUIRE(pairs.size() == 3);
                REQUIRE(containsTupleXOR(pairs, to_persistent_ix(std::make_tuple(a, b))));
                REQUIRE(containsTupleXOR(pairs, to_persistent_ix(std::make_tuple(b, c))));
                REQUIRE(containsTupleXOR(pairs, to_persistent_ix(std::make_tuple(c, a))));
            }

            THEN("Expect 3 unique triples") {
                REQUIRE(triples.size() == 3);
                REQUIRE(containsTupleXOR(triples, to_persistent_ix(std::make_tuple(a, b, c))));
                REQUIRE(containsTupleXOR(triples, to_persistent_ix(std::make_tuple(b, c, a))));
                REQUIRE(containsTupleXOR(triples, to_persistent_ix(std::make_tuple(c, a, b))));
            }

            THEN("Expect 0 unique quadruples") {
                REQUIRE(quadruples.empty());
            }
        }
    }

    testFullyConnected<5>();
    testFullyConnected<7>();
    testFullyConnected<8>();

    GIVEN("Three fully connected graphs, connected by one edge each g1 -- g2 -- g3") {
        auto g1 = fullyConnectedGraph(5);
        auto g2 = fullyConnectedGraph(5);
        auto g3 = fullyConnectedGraph(5);
        auto mapping2 = g1.append(g2, g1.begin(), g2.begin());
        AND_THEN("the graph has the edge [0, mapping(0)]") {
            REQUIRE(g1.containsEdge((g1.begin() + 0).persistent_index(), mapping2.at((g2.begin() + 0).persistent_index().value)));
        }
        auto mapping3 = g1.append(g3, --g1.end(), g3.begin());
        AND_THEN("the graph has the edge [9, mapping(0)]") {
            REQUIRE(g1.containsEdge((g1.begin() + 9).persistent_index(), mapping2.at((g3.begin() + 0).persistent_index().value)));
        }
        REQUIRE(g1.graphDistance(g1.begin()+1, --g1.end()) == 5);
    }

    GIVEN("A fully connected graph g1 of size 5 and a fully connected graph g2 of size 4") {
        auto g1 = fullyConnectedGraph(5);
        auto g1Copy = g1;
        auto g2 = fullyConnectedGraph(4);

        THEN("All graphs are connected") {
            REQUIRE(g1.isConnected());
            REQUIRE(g1Copy.isConnected());
            REQUIRE(g2.isConnected());
        }

        THEN("The distance between all vertices is 1 (unless they are the same vertex)") {
            for(auto it = g1.begin(); it != g1.end(); ++it) {
                for(auto it2 = g1.cbegin_persistent(); it2 != g1.cend_persistent(); ++it2) {
                    if(it.persistent_index() != g1.vertices().persistentIndex(it2)) {
                        REQUIRE(g1.graphDistance(it, it2) == 1);
                    } else {
                        REQUIRE(g1.graphDistance(it2, it) == 0);
                    }
                }
            }
        }

        WHEN("Appending g2 to g1") {
            g1.append(g2);
            THEN("g1 becomes disconnected, while g2 is still connected") {
                REQUIRE_FALSE(g1.isConnected());
                REQUIRE(g2.isConnected());
            }
            THEN("The number of vertices in g1 is 5+4") {
                REQUIRE(g1.vertices().size() == 9);
            }
            THEN("The number of edges is the sum of the number of edges in g1 and g2") {
                REQUIRE(g1.nEdges() == g1Copy.nEdges() + g2.nEdges());
            }
            AND_WHEN("taking the connected components of the joint graph") {
                auto components = g1.connectedComponents();
                THEN("we get back the original g1 and g2") {
                    REQUIRE(components.size() == 2);
                    REQUIRE((components[0].vertices().size() == 5 ^ components[0].vertices().size() == 4) == true);
                    REQUIRE((components[1].vertices().size() == 5 ^ components[1].vertices().size() == 4) == true);
                    if(components[0].vertices().size() == 5) {
                        // components[0] is g1
                        REQUIRE(components[0].nEdges() == g1Copy.nEdges());
                        REQUIRE(components[1].nEdges() == g2.nEdges());
                    } else {
                        // components[0] is g2
                        REQUIRE(components[0].nEdges() == g2.nEdges());
                        REQUIRE(components[1].nEdges() == g1Copy.nEdges());
                    }
                }
            }
        }
        WHEN("Removing a vertex from graph g1") {
            g1.removeVertex(3);
            THEN("The graph is a fully connected graph with 4 vertices") {
                REQUIRE(g1.vertices().size() == 4);
                REQUIRE(g1.vertices().size_persistent() == 5);
                REQUIRE(g1.nEdges() == g2.nEdges());
                for(const auto&[v1, v2] : g1.edges()) {
                    REQUIRE(!g1.vertices().at(v1).deactivated());
                    REQUIRE(!g1.vertices().at(v2).deactivated());
                    REQUIRE(v1 != v2);
                }
            }

            AND_WHEN("Removing a second vertex from g1") {
                g1.removeVertex(2);
                THEN("The graph is a fully connected graph with 3 vertices") {
                    REQUIRE(g1.vertices().size() == 3);
                    REQUIRE(g1.vertices().size_persistent() == 5);
                    REQUIRE(g1.nEdges() == g2.nEdges() - 3);
                }

                AND_WHEN("Removing a third vertex from g1") {
                    g1.removeVertex(1);
                    THEN("The graph is a fully connected graph with 2 vertices") {
                        REQUIRE(g1.vertices().size() == 2);
                        REQUIRE(g1.vertices().size_persistent() == 5);
                        REQUIRE(g1.nEdges() == g2.nEdges() - 3 - 2);
                        for(auto [i1, i2] : g1.edges()) {
                            const auto &v1 = g1.vertices().at(i1);
                            const auto &v2 = g1.vertices().at(i2);
                            REQUIRE(!v1.deactivated());
                            REQUIRE(!v2.deactivated());
                            REQUIRE(std::find(v1.neighbors().begin(), v1.neighbors().end(), i2) != v1.neighbors().end());
                            REQUIRE(std::find(v2.neighbors().begin(), v2.neighbors().end(), i1) != v2.neighbors().end());
                        }
                    }
                }

                AND_WHEN("Appending g2 to g1 forming an edge between 0 in g1 and 3 in g2") {
                    auto mapping = g1.append(g2, 0, 3);
                    THEN("the resulting graph is connected") {
                        REQUIRE(g1.isConnected());
                    }
                    AND_THEN("the graph has the edge [0, mapping(3)]") {
                        REQUIRE(g1.containsEdge((g1.begin() + 0).persistent_index(), mapping.at((g2.begin() + 3).persistent_index().value)));
                    }
                    AND_WHEN("Removing the edge [0, mapping(3)]") {
                        g1.removeEdge((g1.begin() + 0).persistent_index(), mapping.at((g2.begin() + 3).persistent_index().value));
                        THEN("There are two connected components") {
                            REQUIRE(!g1.isConnected());
                            REQUIRE(g1.connectedComponents().size() == 2);
                        }
                        AND_THEN("One component is fully connected with 3 vertices, one with 4 vertices") {
                            auto components = g1.connectedComponents();
                            const auto &gg1 = components.at(0).vertices().size() == 3 ? components.at(0)
                                    : components.at(1);
                            const auto &gg2 = components.at(0).vertices().size() == 3 ? components.at(1)
                                    : components.at(0);
                            REQUIRE(gg1.vertices().size() == 3);
                            REQUIRE(gg1.nEdges() == static_cast<std::size_t>(3 * (3 - 1) / 2));
                            REQUIRE(gg2.vertices().size() == 4);
                            REQUIRE(gg2.nEdges() == g2.nEdges());
                        }
                    }
                }
            }

            AND_WHEN("Appending g2") {
                g1.append(g2);
                THEN("There are two connected components, fully connected graphs with 4 vertices each") {
                    REQUIRE(g1.nEdges() == 2*g2.nEdges());
                    auto components = g1.connectedComponents();
                    REQUIRE(components.size() == 2);
                    REQUIRE(components[0].vertices().size() == 4);
                    REQUIRE(components[1].vertices().size() == 4);
                    REQUIRE(components[0].nEdges() == g2.nEdges());
                    REQUIRE(components[1].nEdges() == g2.nEdges());

                    AND_THEN("The edges list should be consistent with the vertex neighbors") {
                        for(const auto &component : components) {
                            std::size_t nActive = 0;
                            for(std::size_t i = 0; i < component.vertices().size_persistent(); ++i) {
                                auto persistentIx = graphs::PersistentIndex{i};
                                if(!component.vertices().at(persistentIx).deactivated()) {

                                    for(auto neighborIndex : component.vertices().at(i).neighbors()) {
                                        REQUIRE(!component.vertices().at(neighborIndex).deactivated());
                                        auto it = std::find_if(std::begin(component.edges()),
                                                std::end(component.edges()), [persistentIx, neighborIndex](const auto& edge) {
                                                    auto [i1, i2] = edge;
                                                    return (i1 == persistentIx && i2 == neighborIndex) ||
                                                        (i1 == neighborIndex && i2 == persistentIx);
                                        });
                                        REQUIRE(it != component.edges().end());
                                    }

                                    ++nActive;
                                }
                            }
                            REQUIRE(nActive == 4);
                        }
                    }
                }
            }
        }
    }
}
