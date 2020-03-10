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

using PyEdge = std::tuple<PyVertex, PyVertex>;

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

