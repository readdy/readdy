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
 * @file ExportTopologies.cpp
 * @brief Impl file for exporting topology related classes and functionality
 * @author clonker
 * @date 04.02.17
 */

#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <readdy/model/Particle.h>
#include <readdy/model/topologies/GraphTopology.h>
#include <readdy/model/_internal/Util.h>

#include <utility>
#include "PyFunction.h"
#include "../common/ReadableParticle.h"
#include "PyTopology.h"

namespace py = pybind11;

namespace top = readdy::model::top;
namespace reactions = readdy::model::top::reactions;

using rvp = py::return_value_policy;

struct reaction_function_sink {
    std::shared_ptr<py::function> f;
    reaction_function_sink(const py::function& f) : f(std::make_shared<py::function>(f)) {};

    inline reactions::StructuralTopologyReaction::reaction_function::result_type operator()(top::GraphTopology& top) {
        py::gil_scoped_acquire gil;
        PyTopology pyTop (&top);
        auto t = py::cast(&pyTop, py::return_value_policy::automatic_reference);
        auto rv = (*f)(*(t.cast<PyTopology*>()));
        auto pyrecipy = rv.cast<PyRecipe>();
        return pyrecipy.get();
    }
};

struct rate_function_sink {
    std::shared_ptr<py::function> f;
    rate_function_sink(const py::function& f) : f(std::make_shared<py::function>(f)) {};

    inline reactions::StructuralTopologyReaction::rate_function::result_type operator()(const top::GraphTopology& top) {
        py::gil_scoped_acquire gil;
        PyTopology pyTop (&const_cast<top::GraphTopology&>(top));
        auto t = py::cast(&pyTop, py::return_value_policy::automatic_reference);
        auto rv = (*f)(*(t.cast<PyTopology*>()));
        return rv.cast<reactions::StructuralTopologyReaction::rate_function::result_type>();
    }
};


struct NeighborIteratorState {
    PyTopology* top;
    top::Vertex::NeighborList::const_iterator it;
    top::Vertex::NeighborList::const_iterator end;
    bool first_or_done;
};

void exportTopologies(py::module &m) {
    using namespace py::literals;

    py::class_<reaction_function_sink>(m, "ReactionFunction").def(py::init<py::function>());
    py::class_<rate_function_sink>(m, "RateFunction").def(py::init<py::function>());
    py::class_<spatial_rate_function_sink>(m, "SpatialRateFunction").def(py::init<py::function>());

    py::class_<reactions::StructuralTopologyReaction>(m, "StructuralTopologyReaction")
            .def(py::init<std::string, reaction_function_sink, rate_function_sink>())
            .def("rate", &reactions::StructuralTopologyReaction::rate, "topology"_a)
            .def("expects_connected_after_reaction", &reactions::StructuralTopologyReaction::expects_connected_after_reaction)
            .def("expect_connected_after_reaction", &reactions::StructuralTopologyReaction::expect_connected_after_reaction)
            .def("creates_child_topologies_after_reaction", &reactions::StructuralTopologyReaction::creates_child_topologies_after_reaction)
            .def("create_child_topologies_after_reaction", &reactions::StructuralTopologyReaction::create_child_topologies_after_reaction);

    py::class_<PyRecipe>(m, "Recipe")
            .def(py::init<PyTopology>(), R"topdoc(
                 Creates a new reaction recipe.

                 :param topology: The topology for which this recipe should be created.
            )topdoc", "topology"_a)
            .def("change_particle_type", [](PyRecipe &self, const PyVertex &v, const std::string &to) {
                self.get().changeParticleType(v.get(), to);
                return self;
            }, py::return_value_policy::reference_internal)
            .def("change_particle_type", [](PyRecipe &self, std::size_t v, const std::string &to) {
                if (v >= self->topology().graph().nVertices()) {
                    throw std::invalid_argument("vertex index out of bounds (" + std::to_string(v) + " >= "
                                                + std::to_string(self->topology().graph().nVertices()) + ")");
                }
                auto pix = (self->topology().graph().vertices().begin() + v).persistent_index();
                self.get().changeParticleType(pix, to);
                return self;
            }, py::return_value_policy::reference_internal)
            .def("change_particle_position", [](PyRecipe &self, const PyVertex &v, readdy::Vec3 pos) {
                self.get().changeParticlePosition(v.get(), pos);
                return self;
            }, py::return_value_policy::reference_internal)
            .def("change_particle_position", [](PyRecipe &self, std::size_t v, readdy::Vec3 pos) {
                if (v >= self->topology().graph().nVertices()) {
                    throw std::invalid_argument("vertex index out of bounds (" + std::to_string(v) + " >= "
                                                + std::to_string(self->topology().graph().nVertices()) + ")");
                }
                auto pix = (self->topology().graph().vertices().begin() + v).persistent_index();
                self.get().changeParticlePosition(pix, pos);
            })
            .def("add_edge", [](PyRecipe &self, std::size_t v1, std::size_t v2) {
                if (v1 >= self->topology().graph().nVertices() || v2 >= self->topology().graph().nVertices()) {
                    throw std::invalid_argument(
                            fmt::format("At least one vertex index ({}, {}) out of bounds nVertices = {}",
                                        v1, v2, self->topology().graph().nVertices()));
                }
                auto pix1 = (self.get().topology().graph().begin() + v1).persistent_index();
                auto pix2 = (self.get().topology().graph().begin() + v2).persistent_index();
                self.get().addEdge(pix1, pix2);
                return self;
            }, R"topdoc(
                Adds an edge between the given vertices.

                :param v_index1: index of the first vertex, as in `topology.get_graph().get_vertices()`
                :param v_index2: index of the second vertex, as in `topology.get_graph().get_vertices()`
                :return: a reference to this recipe to enable a fluent interface
            )topdoc", "v_index1"_a, "v_index2"_a, py::return_value_policy::reference_internal)
            .def("add_edge", [](PyRecipe &self, const PyVertex &v1, const PyVertex &v2) {
                self.get().addEdge(v1.get(), v2.get());
                return self;
            }, py::return_value_policy::reference_internal)
            .def("append_particle", [](PyRecipe &self, const std::vector<PyVertex> &neighbors,
                    const std::string &type, readdy::Vec3 pos){
                std::vector<top::Graph::PersistentVertexIndex> refs;
                std::transform(neighbors.begin(), neighbors.end(), std::back_inserter(refs), [](const PyVertex& v) {
                    return v.get();
                });
                self.get().appendNewParticle(refs, type, pos);
                return self;
            })
            .def("append_particle", [](PyRecipe &self, const std::vector<std::size_t> &neighbors,
                    const std::string &type, readdy::Vec3 pos){
                std::vector<top::Graph::PersistentVertexIndex> refs;
                std::transform(neighbors.begin(), neighbors.end(), std::back_inserter(refs), [&](auto ix) {
                    return (self.get().topology().graph().begin() + ix).persistent_index();
                });
                self.get().appendNewParticle(refs, type, pos);
                return self;
            })
            .def("remove_edge", [](PyRecipe &self, std::size_t v1, std::size_t v2) -> PyRecipe& {
                if (v1 >= self->topology().graph().nVertices() || v2 >= self->topology().graph().nVertices()) {
                    throw std::invalid_argument(
                            fmt::format("At least one vertex index ({}, {}) out of bounds nVertices = {}",
                                    v1, v2, self->topology().graph().nVertices()));
                }
                auto pix1 = (self.get().topology().graph().begin() + v1).persistent_index();
                auto pix2 = (self.get().topology().graph().begin() + v2).persistent_index();
                self.get().removeEdge(pix1, pix2);
                return self;
            }, R"topdoc(
                Removes an edge between given vertices. Depending on the configuration of the topology reaction, this
                can lead to failed states or multiple sub-topologies.

                :param v_index1: index of the first vertex, as in `topology.get_graph().get_vertices()`
                :param v_index2: index of the second vertex, as in `topology.get_graph().get_vertices()`
                :return: a reference to this recipe to enable a fluent interface
            )topdoc", "v_index1"_a, "v_index2"_a, py::return_value_policy::reference_internal)
            .def("remove_edge", [](PyRecipe &self, const PyVertex &v1, const PyVertex &v2) {
                self->removeEdge(v1.get(), v2.get());
                return self;
            }, py::return_value_policy::reference_internal)
            .def("remove_edge", [](PyRecipe &self, PyEdge edge) {
                self->removeEdge(std::get<0>(edge).get(), std::get<1>(edge).get());
                return self;
            }, R"topdoc(
                Removes an edge between given vertices. Depending on the configuration of the topology reaction, this
                can lead to failed states or multiple sub-topologies.

                :param edge: the edge
                :return: a reference to this recipe to enable a fluent interface
            )topdoc", "edge"_a, py::return_value_policy::reference_internal)
            .def("separate_vertex", [](PyRecipe &self, const std::size_t v) {
                if (v >= self->topology().graph().nVertices()) {
                    throw std::invalid_argument("vertex index out of bounds (" + std::to_string(v) + " >= "
                                                + std::to_string(self->topology().graph().nVertices()) + ")");
                }
                auto pix = (self.get().topology().graph().vertices().begin() + v).persistent_index();
                self.get().separateVertex(pix);
                return self;
            }, R"topdoc(
                Removes all edges from the topology's graph that contain the vertex corresponding to the provided index.

                If no new edge is formed between the given vertex this call, depending on the configuration of the
                reaction, can lead to a failed state or to formation of a topology consisting out of only one particle.
                In the latter case, this call can be followed by a call to `change_particle_type`, where the target
                type is no topology type. Then, no one-particle topology will be formed but the particle will simply
                be emitted and treated as normal particle.

                :param index: The vertex' index with respect to `topology.get_graph().get_vertices()`
                :return: a reference to this recipe to enable a fluent interface
            )topdoc", "index"_a, py::return_value_policy::reference_internal)
            .def("separate_vertex", [](PyRecipe &self, const PyVertex &v) {
                self.get().separateVertex(v.get());
                return self;
            }, py::return_value_policy::reference_internal)
            .def("change_topology_type", [](PyRecipe &self, std::string type) {
                self.get().changeTopologyType(type);
                return self;
                }, R"topdoc(
                Changes the type of the topology to the given type, potentially changing its structural and spatial
                topology reactions.

                :param type: the target type
                :return: a reference to this recipe to enable a fluent interface
            )topdoc","type"_a, py::return_value_policy::reference_internal);

    py::class_<PyTopology>(m, "Topology")
            .def("get_graph", [](PyTopology &self) -> PyGraph { return {&self}; })
            .def_property_readonly("graph", [](PyTopology &self) { return PyGraph{&self}; })
            .def("get_n_particles", [](const PyTopology &self) { return self.get()->graph().nVertices(); })
            .def_property_readonly("n_particles", [](const PyTopology &self) { return self.get()->graph().nVertices(); })
            .def_property_readonly("particles", [](const PyTopology &self) -> std::vector<rpy::ReadableParticle> {
                auto particles = self.get()->fetchParticles();
                std::vector<rpy::ReadableParticle> readableOutput;
                readableOutput.reserve(particles.size());
                std::transform(std::begin(particles), std::end(particles), std::back_inserter(readableOutput),
                               [&self](const auto &particle) { return rpy::ReadableParticle(particle, self.get()->context());});
                return readableOutput;
                }, R"topdoc(
                Retrieves the particles contained in this topology.

                :return: the particles
            )topdoc")
            .def_property_readonly("type", [](const PyTopology &self) {
                return self->context().topologyRegistry().nameOf(self->type());
            })
            .def("particle_type_of_vertex", [](const PyTopology &self, const PyVertex &v) -> std::string {
                return self->context().particleTypes().nameOf(self->particleForVertex(v.get()).type());
            }, R"topdoc(
               Retrieves the particle type corresponding to a vertex.

               :param v: the vertex
               :return: the particle type
            )topdoc")
            .def("position_of_vertex", [](const PyTopology &self, const PyVertex &v) -> readdy::Vec3 {
                return self->particleForVertex(v.get()).pos();
            }, R"topdoc(
                Retrieves the position of the particle corresponding to the given vertex.

                :param v: the vertex
                :return: the position
            )topdoc")
            .def("particle_id_of_vertex", [](const PyTopology &self, const PyVertex &v) -> readdy::ParticleId {
                return self->particleForVertex(v.get()).id();
            }, R"topdoc(
                Retrieves the id of the particle corresponding to the given vertex.

                :param v: the vertex
                :return: the id
            )topdoc")
            .def("configure", [](PyTopology &self) {
                self->configure();
            })
            .def("validate", [](PyTopology &self) {
                self->validate();
            });

    py::class_<PyGraph>(m, "Graph")
            .def("get_vertices", [](PyGraph &self) {
                std::vector<PyVertex> vertices;
                vertices.reserve(self->nVertices());
                for(auto it = self->vertices().begin(); it != self->vertices().end(); ++it) {
                    vertices.emplace_back(self.top(), it.persistent_index());
                }
                return vertices;
            },
            R"topdoc(
                Yields a list of vertices contained in this graph.

                :return: list of vertices
            )topdoc",rvp::copy)
            .def_property_readonly("vertices", [](PyGraph &self) {
                   std::vector<PyVertex> vertices;
                   vertices.reserve(self->nVertices());
                   for(auto it = self->vertices().begin(); it != self->vertices().end(); ++it) {
                       vertices.emplace_back(self.top(), it.persistent_index());
                   }
                   return vertices;
            },
            R"topdoc(
                Yields a list of vertices contained in this graph.

                :return: list of vertices
            )topdoc", rvp::copy)
            .def("get_edges", [](PyGraph &self) -> std::vector<PyEdge> {
                std::vector<PyEdge> edges;
                edges.reserve(self->nEdges());
                for(const auto& edge : self->edges()) {
                    edges.emplace_back(PyVertex(self.top(), std::get<0>(edge)), PyVertex(self.top(), std::get<1>(edge)));
                }
                return edges;
            }, R"topdoc(
                Yields a list of edges contained in this graph.

                :return: list of edges
            )topdoc")
            .def_property_readonly("edges", [](PyGraph &self) {
                std::vector<PyEdge> edges;
                edges.reserve(self->nEdges());
                for(const auto& edge : self->edges()) {
                    edges.emplace_back(PyVertex(self.top(), std::get<0>(edge)), PyVertex(self.top(), std::get<1>(edge)));
                }
                return edges;
            }, R"topdoc(
                Yields a list of edges contained in this graph.

                :return: list of edges
            )topdoc")
            .def("has_edge", [](PyGraph &self, std::size_t v1, std::size_t v2) {
                if(v1 >= self->nVertices() || v2 >= self->nVertices()) {
                    throw std::invalid_argument("vertices out of bounds!");
                }
                auto pix1 = (self->begin() + v1).persistent_index();
                auto pix2 = (self->begin() + v2).persistent_index();
                return self->containsEdge(pix1, pix2);
            }, "vertex_index_1"_a, "vertex_index_2"_a)
            .def("add_edge", [](PyGraph &self, std::size_t v1, std::size_t v2) {
                if(v1 >= self->nVertices() || v2 >= self->nVertices()) {
                    throw std::invalid_argument("vertices out of bounds!");
                }
                auto pix1 = (self->begin() + v1).persistent_index();
                auto pix2 = (self->begin() + v2).persistent_index();
                self.top()->get()->addEdge(pix1, pix2);
            }, "vertex_index_1"_a, "vertex_index_2"_a)
            .def("add_edge", [](PyGraph &self, const PyVertex &v1, const PyVertex &v2) {
                self.top()->get()->addEdge(v1.get(), v2.get());
            });

    py::class_<NeighborIteratorState>(m, "iterator", pybind11::module_local())
            .def("__iter__", [](NeighborIteratorState &s) -> NeighborIteratorState& { return s; })
            .def("__next__", [](NeighborIteratorState &s) -> PyVertex {
                if (!s.first_or_done)
                    ++s.it;
                else
                    s.first_or_done = false;
                if (s.it == s.end) {
                    s.first_or_done = true;
                    throw py::stop_iteration();
                }
                return PyVertex(s.top, *s.it);
            });

    py::class_<PyVertex>(m, "Vertex")
            .def("particle_type", [](const PyVertex &self) {
                auto typeId = self.top()->get()->particleForVertex(self.get()).type();
                return self.top()->get()->context().particleTypes().nameOf(typeId);
            }, R"topdoc(
                Yields this vertex' corresponding particle type.

                :return: the particle type
            )topdoc", rvp::copy)
            .def("neighbors", [](PyVertex &self) {
                const auto& v = self.top()->get()->graph().vertices().at(self.get());
                std::vector<PyVertex> neighbors;
                neighbors.reserve(v.neighbors().size());
                for(auto pix : v.neighbors()) {
                    neighbors.emplace_back(self.top(), pix);
                }
                return neighbors;
            }, R"topdoc(
                Yields this vertex' neighbors.

                :return: this vertex' neighbors.
            )topdoc", rvp::copy)
            .def("__len__", [](const PyVertex &self) {
                const auto& v = self.top()->get()->graph().vertices().at(self.get());
                return v.neighbors().size();
            }, R"topdoc(
                Yields the number of neighbors of this vertex.

                :return: number of neighbors
            )topdoc", rvp::copy)
            .def("__iter__", [](PyVertex &self) {
                return NeighborIteratorState{
                    self.top(),
                    self.top()->get()->graph().vertices().at(self.get()).neighbors().begin(),
                    self.top()->get()->graph().vertices().at(self.get()).neighbors().end(),
                    true
                };
            }, R"topdoc(
                Yields an iterator over this vertex' neighbors.

                :return: the iterator
            )topdoc", py::keep_alive<0, 1>())
            .def_property_readonly("particle_index", [](const PyVertex &v) -> std::size_t{
                return v.top()->get()->graph().vertices().at(v.get())->particleIndex;
            }, R"topdoc(
                Retrieves the particle index for this particle vertex.

                :return: the corresponding particle's index
            )topdoc")
            .def_property("data", [](const PyVertex& v) {
                return v.top()->get()->graph().vertices().at(v.get())->data;
            }, [](PyVertex &v, const std::string& s) {
                v.top()->get()->setVertexData(v.get(), s);
            })
            .def("get", [](PyVertex &v) {
                return PyVertex(v);
            })
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def(py::self < py::self)
            .def(py::self <= py::self)
            .def(py::self > py::self)
            .def(py::self >= py::self)
            .def("__repr__", [](const PyVertex &v) {
                const auto &vertex = v.top()->get()->graph().vertices().at(v.get());
                std::stringstream ss;

                ss << "Vertex[";
                ss << "particleIndex=" << vertex->particleIndex << ", ";
                ss << "data=" << vertex->data << ", ";
                ss << "neighbors=[";
                for(std::size_t i = 0; i < vertex.neighbors().size(); ++i) {
                    if (i > 0) ss << ", ";
                    auto pix = vertex.neighbors().at(i);
                    auto pit = v.top()->get()->graph().vertices().begin_persistent() + pix.value;
                    auto ait = v.top()->get()->graph().vertices().persistent_to_active_iterator(pit);
                    ss << std::distance(v.top()->get()->graph().vertices().begin(), ait);
                }
                ss << "]]";

                return ss.str();
            });
}
