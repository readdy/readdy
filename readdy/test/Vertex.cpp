//
// Created by mho on 11/4/19.
//

#include <catch2/catch_test_macros.hpp>
#include <graphs/Vertex.h>

TEST_CASE("Test Vertex class", "[vertex]") {
    SECTION("With with default ctor") {
        graphs::Vertex vertex (0, 1, 2, 3);
        REQUIRE(vertex.data() == std::make_tuple(0, 1, 2, 3));
        REQUIRE(vertex.neighbors().empty());
    }

    SECTION("With non-default ctor") {
        graphs::Vertex vertex2;
        REQUIRE(vertex2.neighbors().empty());
    }

    SECTION("Vertex with struct") {
        struct A {
            int x;
        };
        A a {
            .x = 20
        };
        graphs::Vertex<A> v (a);
        REQUIRE(v->x == 20);
        v->x = 10;
        REQUIRE(v->x == 10);
    }

    SECTION("Vertices with neighbors") {
        graphs::Vertex<std::size_t> vertex (5);
        REQUIRE(vertex.neighbors().empty());
        
        WHEN("adding a neighbor") {
            vertex.addNeighbor(graphs::PersistentIndex{1});
            THEN("the number of neighbors increases") {
                REQUIRE(vertex.neighbors().size() == 1);
                REQUIRE(vertex.neighbors().front() == graphs::PersistentIndex{1});
            }
            WHEN("adding that neighbor again") {
                vertex.addNeighbor(graphs::PersistentIndex{1});
                THEN("the number of neighbors does not increase") {
                    REQUIRE(vertex.neighbors().size() == 1);
                    REQUIRE(vertex.neighbors().front() == graphs::PersistentIndex{1});
                }
            }

            WHEN("adding a second neighbor") {
                vertex.addNeighbor(graphs::PersistentIndex{2});
                THEN("the number of neighbors is 2") {
                    REQUIRE(vertex.neighbors().size() == 2);
                }
                WHEN("removing the first neighbor") {
                    vertex.removeNeighbor(graphs::PersistentIndex{1});
                    THEN("only the second neighbor remains, now index 1") {
                        REQUIRE(vertex.neighbors().size() == 1);
                        REQUIRE(vertex.neighbors().front() == graphs::PersistentIndex{2});
                    }
                    WHEN("also removing the second neighbor") {
                        vertex.removeNeighbor(graphs::PersistentIndex{2});
                        THEN("no neighbors remain") {
                            REQUIRE(vertex.neighbors().empty());
                        }
                    }
                }
            }
        }
    }

}
