//
// Created by mho on 10/29/19.
//

#pragma once

#include <queue>
#include "../Graph.h"

namespace graphs {

namespace {
template<typename T, size_t... I>
auto reverse_tuple_impl(T t, std::index_sequence<I...>) {
    return std::make_tuple(std::get<sizeof...(I) - 1 - I>(std::forward<T>(t))...);
}

template<typename T>
auto reverse_tuple(T t) {
    return reverse_tuple_impl(std::forward<T>(t), std::make_index_sequence<std::tuple_size<T>::value>());
}
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline Graph<VertexCollection, Vertex, Rest...>::Graph() = default;

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline Graph<VertexCollection, Vertex, Rest...>::Graph(VertexList vertexList) : _vertices(std::move(vertexList)) {
    findEdges([this](const auto& edge) {
        _edges.push_back(edge);
    });
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline Graph<VertexCollection, Vertex, Rest...>::~Graph() = default;

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
template<typename PairCallback, typename TripleCallback, typename QuadrupleCallback>
inline void Graph<VertexCollection, Vertex, Rest...>::findNTuples(const PairCallback &pairCallback,
                                       const TripleCallback &tripleCallback,
                                       const QuadrupleCallback &quadrupleCallback) {
    std::vector<char> visited (_vertices.size_persistent(), false);

    for (std::size_t vertexIndex = 0; vertexIndex < _vertices.size_persistent(); ++vertexIndex) {
        // vertex v1
        auto pvix = PersistentVertexIndex{vertexIndex};
        visited.at(vertexIndex) = true;
        const auto &v1 = vertices().at(PersistentVertexIndex{vertexIndex});
        if(!v1.deactivated()) {
            auto &neighbors = v1.neighbors();
            for (auto neighborIndex : neighbors) {
                // vertex v2 in N(v1)
                if (!visited.at(neighborIndex.value)) {
                    const auto &v2 = _vertices.at(neighborIndex);
                    pairCallback(std::tie(pvix, neighborIndex));
                    for (auto quadIx1 : neighbors) {
                        if (neighborIndex != quadIx1) {
                            // vertex v3 in N(v1)\{v2}
                            for (auto quadIx2 : v2.neighbors()) {
                                if (quadIx2 != pvix && quadIx2 != quadIx1) {
                                    // vertex v4 in N(v2)\{v1, v3}
                                    quadrupleCallback(std::tie(quadIx1, pvix, neighborIndex, quadIx2));
                                }
                            }
                        }
                    }
                }
                for (auto neighborIx2 : neighbors) {
                    if (neighborIx2 != neighborIndex && neighborIx2 < neighborIndex) {
                        tripleCallback(std::tie(neighborIx2, pvix, neighborIndex));
                    }
                }
            }
        }
    }
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline const typename Graph<VertexCollection, Vertex, Rest...>::VertexList &Graph<VertexCollection, Vertex, Rest...>::vertices() const {
    return _vertices;
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::VertexList &Graph<VertexCollection, Vertex, Rest...>::vertices() {
    return _vertices;
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline bool Graph<VertexCollection, Vertex, Rest...>::containsEdge(const Edge &edge) const {
    if(std::find(edges().begin(), edges().end(), edge) != edges().end()) {
        return true;
    } else {
        return std::find(edges().begin(), edges().end(), reverse_tuple(edge)) != edges().end();
    }
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline bool Graph<VertexCollection, Vertex, Rest...>::containsEdge(iterator v1, iterator v2) const {
    return containsEdge(v1.persistent_index(), v2.persistent_index());
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline bool Graph<VertexCollection, Vertex, Rest...>::containsEdge(PersistentVertexIndex v1, PersistentVertexIndex v2) const {
    return containsEdge(std::tie(v1, v2));
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline bool Graph<VertexCollection, Vertex, Rest...>::containsEdge(ActiveVertexIndex v1, ActiveVertexIndex v2) const {
    auto it1 = begin() + v1;
    auto it2 = begin() + v2;
    return containsEdge(it1.persistent_index(), it2.persistent_index());
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::PersistentVertexIndex Graph<VertexCollection, Vertex, Rest...>::addVertex(typename Vertex::data_type data) {
    return _vertices.emplace_back(data);
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::addEdge(iterator it1, iterator it2) {
    addEdge(it1.persistent_index(), it2.persistent_index());
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::addEdge(persistent_iterator it1, persistent_iterator it2) {
    if(it1->deactivated() || it2->deactivated()) {
        throw std::invalid_argument("Tried adding an edge between vertices of which at least one was deactivated.");
    }
    auto ix1 = _vertices.persistentIndex(it1);
    auto ix2 = _vertices.persistentIndex(it2);
    addVertexNeighbor(*it1, ix2);
    addVertexNeighbor(*it2, ix1);
    _edges.push_back(std::make_tuple(ix1, ix2));
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::addEdge(ActiveVertexIndex ix1, ActiveVertexIndex ix2) {
    auto it1 = _vertices.begin() + ix1;
    auto it2 = _vertices.begin() + ix2;
    addVertexNeighbor(*it1, it2.persistent_index());
    addVertexNeighbor(*it2, it1.persistent_index());
    _edges.push_back(std::make_tuple(PersistentVertexIndex{it1.persistent_index()},
                                     PersistentVertexIndex{it2.persistent_index()}));
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::addEdge(PersistentVertexIndex ix1, PersistentVertexIndex ix2) {
    addEdge(_vertices.begin_persistent() + ix1.value, _vertices.begin_persistent() + ix2.value);
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::addEdge(const Edge &edge) {
    addEdge(std::get<0>(edge), std::get<1>(edge));
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
void Graph<VertexCollection, Vertex, Rest...>::removeEdge(iterator it1, iterator it2) {
    removeEdge(it1.persistent_index(), it2.persistent_index());
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
void Graph<VertexCollection, Vertex, Rest...>::removeEdge(persistent_iterator it1, persistent_iterator it2) {
    if(it1->deactivated() || it2->deactivated()) {
        throw std::invalid_argument("Tried removing an edge between vertices of which at least one was deactivated.");
    }
    auto ix1 = _vertices.persistentIndex(it1);
    auto ix2 = _vertices.persistentIndex(it2);
    auto it = std::find_if(_edges.begin(), _edges.end(), [ix1, ix2](const auto& edge) {
        const auto &[e1, e2] = edge;
        return (e1 == ix1 && e2 == ix2) || (e1 == ix2 && e2 == ix1);
    });
    if(it != edges().end()) {
        removeVertexNeighbor(*it1, ix2);
        removeVertexNeighbor(*it2, ix1);
        _edges.erase(it);
    }/* else {
        throw std::invalid_argument("Tried to remove non-existing edge!");
    }*/
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::removeEdge(ActiveVertexIndex ix1, ActiveVertexIndex ix2) {
    removeEdge(begin() + ix1, begin() + ix2);
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::removeEdge(PersistentVertexIndex ix1, PersistentVertexIndex ix2) {
    removeEdge(_vertices.begin_persistent() + ix1.value, _vertices.begin_persistent() + ix2.value);
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::removeEdge(const Edge &edge) {
    removeEdge(std::get<0>(edge), std::get<1>(edge));
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
void Graph<VertexCollection, Vertex, Rest...>::removeVertex(iterator it) {
    removeVertex(it.to_persistent());
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
void Graph<VertexCollection, Vertex, Rest...>::removeVertex(persistent_iterator it) {
    auto ix = _vertices.persistentIndex(it);
    removeNeighborsEdges(ix);
    _vertices.erase(it);
    _edges.erase(std::remove_if(_edges.begin(), _edges.end(), [ix](const auto &edge) {
        return std::get<0>(edge) == ix || std::get<1>(edge) == ix;
    }), _edges.end());
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::removeVertex(PersistentVertexIndex ix) {
    removeVertex(_vertices.begin_persistent() + ix.value);
}


template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::removeVertex(ActiveVertexIndex ix) {
    auto it = begin() + ix;
    removeVertex(it.persistent_index());
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline const std::vector<typename Graph<VertexCollection, Vertex, Rest...>::Edge> &Graph<VertexCollection, Vertex, Rest...>::edges() const {
    return _edges;
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
template<typename PairCallback>
inline void Graph<VertexCollection, Vertex, Rest...>::findEdges(const PairCallback &edgeCallback) const {
    for(auto it = _vertices.begin(); it != _vertices.end(); ++it) {
        const auto& ix = it.persistent_index();
        for(auto neighbor : it->neighbors()) {
            if(neighbor < it.persistent_index()) {
                const Edge e{std::tie(ix, neighbor)};
                edgeCallback(e);
            }
        }
    }
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline std::tuple<std::vector<typename Graph<VertexCollection, Vertex, Rest...>::Edge>,
                  std::vector<typename Graph<VertexCollection, Vertex, Rest...>::Path3>,
                  std::vector<typename Graph<VertexCollection, Vertex, Rest...>::Path4>> Graph<VertexCollection, Vertex, Rest...>::findNTuples() {
    auto tuple = std::make_tuple(std::vector<Edge>(), std::vector<Path3>(), std::vector<Path4>());
    findNTuples([&](const Edge &edge) {
        std::get<0>(tuple).push_back(edge);
    }, [&](const Path3 &path3) {
        std::get<1>(tuple).push_back(path3);
    }, [&](const Path4 &path4) {
        std::get<2>(tuple).push_back(path4);
    });
    return tuple;
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::addVertexNeighbor(Vertex &v1, PersistentVertexIndex v2) {
    v1.addNeighbor(v2);
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::removeVertexNeighbor(Vertex &v1, PersistentVertexIndex v2) {
    v1.removeNeighbor(v2);
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
template<typename T1, typename T2>
std::int32_t Graph<VertexCollection, Vertex, Rest...>::graphDistance(T1 it1, T2 it2) const {
    const_persistent_iterator it1Persistent = toPersistentIterator(it1);
    const_persistent_iterator it2Persistent = toPersistentIterator(it2);
    PersistentIndex ixSource = _vertices.persistentIndex(it1Persistent);
    PersistentIndex ixTarget = _vertices.persistentIndex(it2Persistent);

    std::vector<char> visited (_vertices.size_persistent(), false);
    std::vector<std::int32_t> dists (_vertices.size_persistent(), 0);

    std::queue<PersistentVertexIndex> unvisited;
    unvisited.push(ixSource);
    visited[ixSource.value] = true;

    while(!unvisited.empty()) {
        auto ix = unvisited.front();
        unvisited.pop();
        auto dCurr = dists[ix.value];

        const auto &vertex = _vertices.at(ix);
        for(auto neighbor : vertex.neighbors()) {
            if(visited[neighbor.value]) {
                continue;
            }

            dists[neighbor.value] = dCurr + 1;
            unvisited.push(neighbor);
            visited[neighbor.value] = true;

            if(neighbor == ixTarget) {
                break;
            }
        }
    }

    if(visited[ixTarget.value]) {
        return dists[ixTarget.value];
    }
    return -1;
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline bool Graph<VertexCollection, Vertex, Rest...>::isConnected() const {
    if(_vertices.empty()) return true;

    std::vector<char> visited (_vertices.size_persistent(), false);

    std::vector<PersistentVertexIndex> unvisited;
    unvisited.reserve(_vertices.size_persistent());
    unvisited.emplace_back([this](){
        auto it = _vertices.begin();
        if(it != _vertices.end()) {
            return it.persistent_index();
        }
        throw std::logic_error("vertices is not empty but could not find an active vertex");
    }());

    std::size_t nVisited = 0;
    while(!unvisited.empty()) {
        auto vertexIndex = unvisited.back();
        const auto &vertex = _vertices.at(vertexIndex);
        unvisited.pop_back();
        if(!visited.at(vertexIndex.value)) {
            visited.at(vertexIndex.value) = true;
            ++nVisited;
            for (auto neighbor : vertex.neighbors()) {
                if (!visited.at(neighbor.value)) {
                    unvisited.emplace_back(neighbor);
                }
            }
        }
    }
    return nVisited == _vertices.size();
}

template<template <typename...>class VertexCollection, typename Vertex, typename... Rest>
inline auto Graph<VertexCollection, Vertex, Rest...>::connectedComponents() -> std::vector<Graph>{
    std::vector<VertexList> subVertexLists {};
    {
        std::vector<std::vector<PersistentVertexIndex>> components {};
        std::vector<char> visited (_vertices.size_persistent(), false);

        for(std::size_t ix = 0; ix < _vertices.size_persistent(); ++ix) {
            if(!_vertices.at(PersistentVertexIndex{ix}).deactivated() && !visited.at(ix)) {
                // got a new component
                components.emplace_back();
                subVertexLists.emplace_back();

                auto& component = components.back();

                std::vector<PersistentVertexIndex> unvisitedInComponent;
                unvisitedInComponent.emplace_back(PersistentVertexIndex{ix});
                while (!unvisitedInComponent.empty()) {
                    auto vertexIndex = unvisitedInComponent.back();
                    unvisitedInComponent.pop_back();
                    if (!visited.at(vertexIndex.value)) {
                        visited.at(vertexIndex.value) = true;
                        component.emplace_back(vertexIndex);
                        for (auto neighbor : _vertices.at(vertexIndex).neighbors()) {
                            if (!visited.at(neighbor.value)) {
                                unvisitedInComponent.emplace_back(neighbor);
                            }
                        }
                    }
                }
            }
        }

        {
            // transfer vertices
            auto itComponents = components.begin();
            auto itSubLists = subVertexLists.begin();

            // mapping (previous vertex index) -> (new vertex index)
            std::vector<std::vector<PersistentVertexIndex>> reverseMappings;
            {
                reverseMappings.reserve(components.size());
                for(const auto &component : components) {
                    std::vector<PersistentVertexIndex> reverseMapping {};
                    reverseMapping.resize(_vertices.size_persistent());
                    for(std::size_t i = 0; i < component.size(); ++i) {
                        reverseMapping[component.at(i).value] = PersistentVertexIndex{i};
                    }
                    reverseMappings.push_back(std::move(reverseMapping));
                }
            }

            auto itReverseMappings = reverseMappings.begin();

            for(; itComponents != components.end(); ++itComponents, ++itSubLists, ++itReverseMappings) {
                for(std::size_t i = 0; i < (*itComponents).size(); ++i) {
                    auto previousVertexIndex = itComponents->at(i);
                    itSubLists->push_back(_vertices.at(previousVertexIndex));
                    for(auto &neighborIndex : (itSubLists->end()-1)->neighbors()) {
                        neighborIndex = itReverseMappings->at(neighborIndex.value);
                    }
                }
            }
        }
    }

    std::vector<Graph<VertexCollection, Vertex, Rest...>> subGraphs {};
    subGraphs.reserve(subVertexLists.size());
    {
        for (auto &subVertexList : subVertexLists) {
            subGraphs.emplace_back(std::move(subVertexList));
        }
    }
    return std::move(subGraphs);
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline void Graph<VertexCollection, Vertex, Rest...>::removeNeighborsEdges(PersistentVertexIndex ix) {
    auto &vertex = *(_vertices.begin_persistent() + ix.value);
    std::for_each(std::begin(vertex.neighbors()), std::end(vertex.neighbors()), [this, ix](const auto neighbor) {
        removeVertexNeighbor(_vertices.at(neighbor), ix);
    });
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline std::size_t Graph<VertexCollection, Vertex, Rest...>::nEdges() const {
    return edges().size();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline std::vector<typename Graph<VertexCollection, Vertex, Rest...>::PersistentVertexIndex> Graph<VertexCollection, Vertex, Rest...>::append(const Graph &other) {
    std::vector<PersistentVertexIndex> indexMapping;
    indexMapping.resize(other.vertices().size_persistent());
    // insert vertices
    {
        std::size_t ix{0};
        for(auto it = other.vertices().begin_persistent(); it != other.vertices().end_persistent(); ++it, ++ix) {
            if (!it->deactivated()) {
                auto persistentIndex = addVertex(it->data());
                indexMapping.at(ix) = persistentIndex;
            }
        }
    }
    // insert edges
    {
        std::size_t ix{0};
        for(auto it = other.vertices().begin_persistent(); it != other.vertices().end_persistent(); ++it, ++ix) {
            if(!it->deactivated()) {
                auto newIndex = indexMapping.at(ix);
                for(auto neighborIndex : it->neighbors()) {
                    // only add once since neighbors are bidirectional
                    if(ix < neighborIndex.value) {
                        addEdge(newIndex, indexMapping.at(neighborIndex.value));
                    }
                }
            }
        }
    }
    return std::move(indexMapping);
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline std::vector<typename Graph<VertexCollection, Vertex, Rest...>::PersistentVertexIndex> Graph<VertexCollection, Vertex, Rest...>::append(
        const Graph &other,
        ActiveVertexIndex edgeIndexThis,
        ActiveVertexIndex edgeIndexOther) {
    return append(other, (begin() + edgeIndexThis).persistent_index(), (other.begin() + edgeIndexOther).persistent_index());
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline std::vector<typename Graph<VertexCollection, Vertex, Rest...>::PersistentVertexIndex> Graph<VertexCollection, Vertex, Rest...>::append(
        const Graph &other,
        iterator itThis,
        iterator itOther) {
    return append(other, itThis.persistent_index(), itOther.persistent_index());
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline std::vector<typename Graph<VertexCollection, Vertex, Rest...>::PersistentVertexIndex> Graph<VertexCollection, Vertex, Rest...>::append(
        const Graph &other,
        persistent_iterator itThis,
        persistent_iterator itOther) {
    return append(other, _vertices.persistentIndex(itThis), _vertices.persistentIndex(itOther));
}


template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline std::vector<typename Graph<VertexCollection, Vertex, Rest...>::PersistentVertexIndex> Graph<VertexCollection, Vertex, Rest...>::append(
        const Graph &other,
        PersistentVertexIndex edgeIndexThis,
        PersistentVertexIndex edgeIndexOther) {
    auto mapping = append(other);
    addEdge(edgeIndexThis, mapping.at(edgeIndexOther.value));
    return std::move(mapping);
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::VertexList::size_type Graph<VertexCollection, Vertex, Rest...>::nVertices() const {
    return _vertices.size();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline std::string Graph<VertexCollection, Vertex, Rest...>::gexf() const {
    std::ostringstream ss;
    ss << R"(<?xml version="1.0" encoding="UTF-8"?>)";
    ss << R"(<gexf xmlns="http://www.gexf.net/1.2draft" version="1.2">)";
    ss << R"(<graph mode="static" defaultedgetype="undirected">)";
    {
        ss << "<nodes>";
        std::size_t id = 0;
        for (const auto &v : _vertices) {
            if(!v.deactivated()) {
                ss << fmt::format(R"(<node id="{}" />)", id);
            }
            ++id;
        }
        ss << "</nodes>";
    }
    {
        ss << "<edges>";
        std::size_t id = 0;
        for(auto [i1, i2] : _edges) {
            ss << fmt::format(R"(<edge id="{}" source="{}" target="{}" />)", id, i1.value, i2.value);
            ++id;
        }
        ss << "</edges>";
    }
    ss << "</graph>";
    ss << "</gexf>";
    return ss.str();
}


template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::persistent_iterator Graph<VertexCollection, Vertex, Rest...>::begin_persistent() {
    return _vertices.begin_persistent();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::const_persistent_iterator Graph<VertexCollection, Vertex, Rest...>::begin_persistent() const {
    return _vertices.cbegin_persistent();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::const_persistent_iterator Graph<VertexCollection, Vertex, Rest...>::cbegin_persistent() const {
    return _vertices.cbegin_persistent();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::persistent_iterator Graph<VertexCollection, Vertex, Rest...>::end_persistent() {
    return _vertices.end_persistent();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::const_persistent_iterator Graph<VertexCollection, Vertex, Rest...>::end_persistent() const {
    return _vertices.cend_persistent();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::const_persistent_iterator Graph<VertexCollection, Vertex, Rest...>::cend_persistent() const {
    return _vertices.cend_persistent();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::VertexList::iterator Graph<VertexCollection, Vertex, Rest...>::begin() {
    return _vertices.begin();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::VertexList::const_iterator Graph<VertexCollection, Vertex, Rest...>::begin() const {
    return _vertices.cbegin();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::VertexList::const_iterator Graph<VertexCollection, Vertex, Rest...>::cbegin() const {
    return _vertices.cbegin();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::VertexList::iterator Graph<VertexCollection, Vertex, Rest...>::end() {
    return _vertices.end();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::VertexList::const_iterator Graph<VertexCollection, Vertex, Rest...>::end() const {
    return _vertices.cend();
}

template<template<typename...> class VertexCollection, typename Vertex, typename... Rest>
inline typename Graph<VertexCollection, Vertex, Rest...>::VertexList::const_iterator Graph<VertexCollection, Vertex, Rest...>::cend() const {
    return _vertices.cend();
}

}
