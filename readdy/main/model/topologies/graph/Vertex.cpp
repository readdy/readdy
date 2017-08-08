/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
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
 *
 *
 * @file Vertex.cpp
 * @brief 
 * @author clonker
 * @date 4/14/17
 */

#include <readdy/model/topologies/graph/Vertex.h>
#include <readdy/model/topologies/graph/Graph.h>

namespace readdy {
namespace model {
namespace top {
namespace graph {


VertexRef::VertexRef(Vertex::vertex_ptr it) : it(it), graph(nullptr) {}

VertexRef::VertexRef(Graph *const graph, const Vertex::label_t &label) : graph(graph), label(label) {
    if(graph->vertexLabelMapping().find(label) == graph->vertexLabelMapping().end()) {
        throw std::invalid_argument("tried creating a vertex ref with label \"" +
                                            label + "\" which was not contained in the graph vertex label mapping");
    }
}

Vertex &VertexRef::operator*() {
    if(graph == nullptr) {
        return it.operator*();
    }
    return graph->namedVertex(label);
}

Vertex *VertexRef::operator->() {
    if(graph == nullptr) {
        return it.operator->();
    }
    return graph->vertexLabelMapping().at(label).operator->();
}

Vertex::vertex_ptr& VertexRef::data() {
    if(graph == nullptr) {
        return it;
    }
    return graph->vertexLabelMapping().at(label).data();
}

bool VertexRef::operator==(const VertexRef &rhs) const {
    return data() == rhs.data();
}

const Vertex::vertex_ptr &VertexRef::data() const {
    if(graph == nullptr) {
        return it;
    }
    return graph->vertexLabelMapping().at(label).data();
}

bool VertexRef::operator!=(const VertexRef &rhs) const {
    return !(rhs == *this);
}

const Vertex *VertexRef::operator->() const {
    if(graph == nullptr) {
        return it.operator->();
    }
    return graph->vertexLabelMapping().at(label).operator->();
}

const Vertex &VertexRef::operator*() const {
    if(graph == nullptr) {
        return it.operator*();
    }
    return graph->namedVertex(label);
}

std::ostream &operator<<(std::ostream &os, const VertexRef &vertex) {
    if(vertex.graph != nullptr) {
         os << "Vertex(label=" << vertex.label <<")";
    } else {
        os << *vertex;
    }
    return os;
}

bool VertexCRef::operator==(const VertexCRef &rhs) const {
    return data() == rhs.data();
}

bool VertexCRef::operator!=(const VertexCRef &rhs) const {
    return !(rhs == *this);
}

Vertex::vertex_cptr VertexCRef::data() const {
    if(graph != nullptr) {
        return graph->vertexLabelMapping().at(label).data();
    }
    return it;
}

const Vertex &VertexCRef::operator*() const {
    if(graph == nullptr) {
        return it.operator*();
    }
    return graph->namedVertex(label);
}

const Vertex *VertexCRef::operator->() const {
    if(graph == nullptr) {
        return it.operator->();
    }
    const auto& ref = graph->vertexLabelMapping().at(label);
    return ref.operator->();
}


VertexCRef::VertexCRef(Vertex::vertex_cptr it) : it(it), graph(nullptr){

}

VertexCRef::VertexCRef(Vertex::vertex_ptr it) : it(it), graph(nullptr) {}

VertexCRef::VertexCRef(const Graph *const graph, const Vertex::label_t &label) : graph(graph), label(label) {

}

VertexCRef::VertexCRef(const VertexRef &ref) : VertexCRef(ref.data()) {}

std::ostream &operator<<(std::ostream &os, const VertexCRef &vertex) {
    if(vertex.graph != nullptr) {
        os << "Vertex(label=" << vertex.label <<")";
    } else {
        os << *vertex;
    }
    return os;
}

VertexRef::VertexRef() = default;

const Vertex::label_t &Vertex::label() const {
    return _label;
}

Vertex::label_t &Vertex::label() {
    return _label;
}
}
}
}
}