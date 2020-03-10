//
// Created by mho on 3/10/20.
//

#pragma once


class PyTopology {
        public:
        explicit PyTopology(readdy::model::top::GraphTopology* topology) : topology(topology) {}

        readdy::model::top::GraphTopology* get() {
            return topology;
        }

        [[nodiscard]] const readdy::model::top::GraphTopology* get() const {
            return topology;
        }

        auto *operator->() {return topology;}
        const auto*operator->() const {return topology;}
        private:
        readdy::model::top::GraphTopology* topology;
};

class PyGraph {
        public:
        PyGraph(PyTopology* parent) : parent(parent) {}

        [[nodiscard]] const readdy::model::top::Graph* get() const {
            return &parent->get()->graph();
        }

        const auto *operator->() const { return get(); }

        PyTopology* top() { return parent; }
        private:
        PyTopology* parent;
};

class PyVertex {
        public:
        using VertexRef = readdy::model::top::Graph::PersistentVertexIndex;
        PyVertex() : parent(nullptr), vertex({0}) {}
        PyVertex(PyTopology* parent, VertexRef vertex) : parent(parent), vertex(vertex) {}

        [[nodiscard]] VertexRef get() const { return vertex; }
        PyTopology* top() { return parent; }

        const PyTopology* top() const {return parent;}
        private:
        PyTopology* parent;
        VertexRef vertex;
};

class PyEdge {
        public:
        using EdgeRef = readdy::model::top::Graph::Edge;

        PyEdge() : vertices(std::make_tuple(PyVertex(), PyVertex())) {}

        PyEdge(PyVertex v1, PyVertex v2) : vertices(std::make_tuple(v1, v2)) {}

        PyEdge(std::tuple<PyVertex, PyVertex> edge) : vertices(std::move(edge)) {}

        EdgeRef get() {
            return std::make_tuple(std::get<0>(vertices).get(), std::get<1>(vertices).get());
        }

        private:
        std::tuple<PyVertex, PyVertex> vertices;
};

class PyRecipe {
        public:
        PyRecipe(PyTopology top) : top(top), recipe(std::make_shared<readdy::model::top::reactions::Recipe>(*top.get())) {}

        [[nodiscard]] readdy::model::top::reactions::Recipe& get() { return *recipe; }

        readdy::model::top::reactions::Recipe* operator->() {
            return recipe.get();
        }

        private:
        PyTopology top;
        std::shared_ptr<readdy::model::top::reactions::Recipe> recipe;
};

