//
// Created by mho on 3/10/20.
//

#pragma once

#include <pybind11/stl.h>
#include <readdy/model/topologies/GraphTopology.h>

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

        bool operator==(const PyVertex &other) const {
            if (parent != nullptr && other.parent != nullptr) {
                return parent == other.parent && vertex == other.vertex;
            }
            return false;
        }

        bool operator!=(const PyVertex &other) const {
            if (parent != nullptr && other.parent != nullptr) {
                return !operator==(other);
            }
            return false;
        }

        bool operator<(const PyVertex &other) const {
            if (parent != nullptr && other.parent != nullptr) {
                return parent == other.parent && vertex < other.vertex;
            }
            return false;
        }

        bool operator>=(const PyVertex &other) const {
            if (parent != nullptr && other.parent != nullptr) {
                return !(*this < other);
            }
            return false;
        }

        bool operator>(const PyVertex &other) const {
            if (parent != nullptr && other.parent != nullptr) {
                return *this >= other && *this != other;
            }
            return false;
        }

        bool operator<=(const PyVertex &other) const {
            if (parent != nullptr && other.parent != nullptr) {
                return !(*this > other);
            }
            return false;
        }

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

struct spatial_rate_function_sink {
    std::shared_ptr<pybind11::function> f;
    explicit spatial_rate_function_sink(const pybind11::function& f) : f(std::make_shared<pybind11::function>(f)) {};

    inline readdy::model::top::reactions::SpatialTopologyReaction::rate_function::result_type operator()(
            const readdy::model::top::GraphTopology& top1,
            const readdy::model::top::GraphTopology& top2
    ) {
        pybind11::gil_scoped_acquire gil;
        PyTopology pyTop1 (&const_cast<readdy::model::top::GraphTopology&>(top1));
        PyTopology pyTop2 (&const_cast<readdy::model::top::GraphTopology&>(top2));
        auto t1 = pybind11::cast(&pyTop1, pybind11::return_value_policy::automatic_reference);
        auto t2 = pybind11::cast(&pyTop2, pybind11::return_value_policy::automatic_reference);
        auto rv = (*f)(*(t1.cast<PyTopology*>()), *(t2.cast<PyTopology*>()));
        return rv.cast<readdy::model::top::reactions::SpatialTopologyReaction::rate_function::result_type>();
    }
};