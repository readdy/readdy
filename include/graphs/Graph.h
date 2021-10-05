//
// Created by mho on 10/28/19.
//

#pragma once

#include <list>
#include <algorithm>
#include <vector>
#include <sstream>

#include <fmt/format.h>

#include "IndexPersistentVector.h"

namespace graphs {

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
class Graph {
public:

    using VertexList = VertexCollection<Vertex, Rest...>;
    // this points a contiguous list with automatic jumps over the blanks, i.e., behaves like an ordinary vector
    // the index invalidates as soon as vertices are removed, therefore should not be stored but only used
    // instantaneously
    using ActiveVertexIndex = typename VertexList::size_type;
    // persistent index, never invalidates (unless pointed-to vertex is removed)
    using PersistentVertexIndex = typename VertexList::persistent_index_t;

    using Edge = std::tuple<PersistentVertexIndex, PersistentVertexIndex>;
    using Path2 = Edge;
    using Path3 = std::tuple<PersistentVertexIndex, PersistentVertexIndex, PersistentVertexIndex>;
    using Path4 = std::tuple<PersistentVertexIndex, PersistentVertexIndex, PersistentVertexIndex, PersistentVertexIndex>;

    using iterator = typename VertexList::iterator;
    using const_iterator = typename VertexList::const_iterator;

    using persistent_iterator = typename VertexList::persistent_iterator;
    using const_persistent_iterator = typename VertexList::const_persistent_iterator;

    Graph();

    explicit Graph(VertexList vertexList);

    virtual ~Graph();

    Graph(const Graph &) = default;

    Graph &operator=(const Graph &) = default;

    Graph(Graph &&) = default;

    Graph &operator=(Graph &&) = default;

    iterator begin();
    const_iterator begin() const;
    const_iterator cbegin() const;

    iterator end();
    const_iterator end() const;
    const_iterator cend() const;

    persistent_iterator begin_persistent();
    const_persistent_iterator begin_persistent() const;
    const_persistent_iterator cbegin_persistent() const;

    persistent_iterator end_persistent();
    const_persistent_iterator end_persistent() const;
    const_persistent_iterator cend_persistent() const;

    const VertexList &vertices() const;

    VertexList &vertices();

    typename VertexList::size_type nVertices() const;

    bool containsEdge(const Edge &edge) const;

    bool containsEdge(PersistentVertexIndex v1, PersistentVertexIndex v2) const;

    bool containsEdge(ActiveVertexIndex v1, ActiveVertexIndex v2) const;

    bool containsEdge(iterator it1, iterator it2) const;

    bool containsEdge(persistent_iterator it1, persistent_iterator it2) const;

    PersistentVertexIndex addVertex(typename Vertex::data_type data = {});

    void addEdge(iterator it1, iterator it2);

    void addEdge(persistent_iterator it1, persistent_iterator it2);

    void addEdge(PersistentVertexIndex ix1, PersistentVertexIndex ix2);

    void addEdge(ActiveVertexIndex ix1, ActiveVertexIndex ix2);

    void addEdge(const Edge &edge);

    void removeEdge(iterator it1, iterator it2);

    void removeEdge(persistent_iterator it1, persistent_iterator it2);

    void removeEdge(ActiveVertexIndex ix1, ActiveVertexIndex ix2);

    void removeEdge(PersistentVertexIndex ix1, PersistentVertexIndex ix2);

    void removeEdge(const Edge &edge);

    void removeVertex(iterator it);

    void removeVertex(persistent_iterator it);

    void removeVertex(ActiveVertexIndex ix);

    void removeVertex(PersistentVertexIndex ix);

    bool isConnected() const;

    /**
     * Find shortest distance between two vertices in a graph.
     *
     * @tparam T1 vertex1 type
     * @tparam T2 vertex2 type
     * @param it1 vertex1
     * @param it2 vertex2
     * @param leq early stopping, if -1 ignore     * @return shortest distance or -1 if there is no path
     */
    template<typename T1, typename T2>
    std::int32_t graphDistance(T1 it1, T2 it2) const;

    const std::vector<Edge> &edges() const;

    std::size_t nEdges() const;

    template<typename TupleCallback>
    void findEdges(const TupleCallback &edgeCallback) const;

    template<typename TupleCallback, typename TripleCallback, typename QuadrupleCallback>
    void findNTuples(const TupleCallback &pairCallback,
                     const TripleCallback &tripleCallback,
                     const QuadrupleCallback &quadrupleCallback);

    std::tuple<std::vector<Edge>, std::vector<Path3>, std::vector<Path4>> findNTuples();

    /**
     * Returns the connected components in terms of a list of new graph objects
     * @return connected components
     */
    std::vector<Graph> connectedComponents();

    /**
     * Appends the graph `other` to this graph. No edge is introduced, this graph will have at least two connected
     * components afterwards. Returns an index mapping for the `other` graph which (for the active vertices) contains
     * the index of the vertex in this graph, i.e., `newIndex = mapping[oldIndex]`.
     * @param other the other graph
     * @return index mapping
     */
    std::vector<PersistentVertexIndex> append(const Graph &other);

    /**
     * Appends other graph to this one and introduces an edge connecting the two.
     * @param other other graph
     * @param edgeIndexThis index pointing to a live vertex in this graph
     * @param edgeIndexOther index pointing to a live vertex in the other graph
     * @return index mapping
     */
    std::vector<PersistentVertexIndex> append(const Graph &other, iterator itThis, iterator itOther);
    std::vector<PersistentVertexIndex> append(const Graph &other, persistent_iterator itThis, persistent_iterator itOther);
    std::vector<PersistentVertexIndex> append(const Graph &other, PersistentVertexIndex edgeIndexThis, PersistentVertexIndex edgeIndexOther);
    std::vector<PersistentVertexIndex> append(const Graph &other, ActiveVertexIndex edgeIndexThis, ActiveVertexIndex edgeIndexOther);

    std::string gexf() const;

private:
    VertexList _vertices{};
    std::vector<Edge> _edges {};

    void removeNeighborsEdges(PersistentVertexIndex ix);

    /**
     * this has always to be called for both v1 and v2 (symmetric neighborship)
     * @tparam debug
     * @param v1
     * @param v2
     */
    void addVertexNeighbor(Vertex& v1, PersistentVertexIndex v2);

    /**
     * this has always to be called for both v1 and v2 (symmetric neighborship)
     * @tparam debug
     * @param v1
     * @param v2
     */
    void removeVertexNeighbor(Vertex& v1, PersistentVertexIndex v2);

    template<typename T>
    auto toPersistentIterator(T it) const {
        static_assert(std::is_same_v<T, const_persistent_iterator> || std::is_same_v<T, const_iterator> ||
                std::is_same_v<T, persistent_iterator> || std::is_same_v<T, iterator> ||
                std::is_same_v<T, PersistentVertexIndex> || std::is_same_v<T, ActiveVertexIndex>);
        if constexpr (std::is_same_v<T, const_persistent_iterator>) {
            return it;
        } else if constexpr (std::is_same_v<T, const_iterator>) {
            return it.to_persistent();
        } else if constexpr (std::is_same_v<T, persistent_iterator>) {
            return it;
        } else if constexpr (std::is_same_v<T, iterator>) {
            return it.to_persistent();
        } else if constexpr (std::is_same_v<T, PersistentVertexIndex>) {
            return _vertices.begin_persistent() + it.value;
        } else if constexpr (std::is_same_v<T, ActiveVertexIndex>) {
            return (begin() + it).to_persistent();
        }
    }

};

}

#include "bits/Graph_detail.h"
