# coding=utf-8

# Copyright © 2018 Computational Molecular Biology Group,
#                  Freie Universität Berlin (GER)
#
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the
# following conditions are met:
#  1. Redistributions of source code must retain the above
#     copyright notice, this list of conditions and the
#     following disclaimer.
#  2. Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials
#     provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of
#     its contributors may be used to endorse or promote products
#     derived from this software without specific
#     prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
Created on 02.09.18

@author: clonker
"""

from .._internal.readdybinding.api.top import Recipe as _Recipe
from readdy.api.utils import vec3_of as _v3_of

class StructuralReactionRecipe(object):

    def __init__(self, topology):
        """
        Creates a new reaction recipe.

        :param topology: The topology for which this recipe should be created.
        """
        self._recipe = _Recipe(topology)

    def change_particle_type(self, vertex, type_to):
        """
        Changes the particle type of the to vertex associated particle to the given type.
        :param vertex: the vertex or its index (obtainable from `topology.get_graph().get_vertices()`)
        :param type_to: the target particle type
        :return: a reference to this recipe to enable a fluent interface
        """
        self._recipe.change_particle_type(vertex, type_to)
        return self

    def change_particle_position(self, vertex, new_position):
        """
        Changes the particle position of the to vertex associated particle to the given position.
        :param vertex: the vertex or its index (obtainable from `topology.get_graph().get_vertices()`)
        :param type_to: the target particle type
        :return: a reference to this recipe to enable a fluent interface
        """
        self._recipe.change_particle_position(vertex, _v3_of(new_position))
        return self

    def add_edge(self, vertex1, vertex2):
        """
        Adds an edge between the given vertices.

        :param vertex1: vertex or index of the first vertex, as in `topology.get_graph().get_vertices()`
        :param vertex2: vertex or index of the second vertex, as in `topology.get_graph().get_vertices()`
        :return: a reference to this recipe to enable a fluent interface
        """
        self._recipe.add_edge(vertex1, vertex2)
        return self

    def remove_edge(self, vertex1, vertex2):
        """
        Removes an edge between given vertices. Depending on the configuration of the topology reaction, this
        can lead to failed states or multiple sub-topologies.

        :param vertex1: vertex or index of the first vertex, as in `topology.get_graph().get_vertices()`
        :param vertex2: vertex or index of the second vertex, as in `topology.get_graph().get_vertices()`
        :return: a reference to this recipe to enable a fluent interface
        """
        self._recipe.remove_edge(vertex1, vertex2)
        return self

    def separate_vertex(self, vertex):
        """
        Removes all edges from the topology's graph that contain the vertex corresponding to the provided index.

        If no new edge is formed between the given vertex this call, depending on the configuration of the
        reaction, can lead to a failed state or to formation of a topology consisting out of only one particle.
        In the latter case, this call can be followed by a call to `change_particle_type`, where the target
        type is no topology type. Then, no one-particle topology will be formed but the particle will simply
        be emitted and treated as normal particle.

        :param index: vertex or the vertex' index with respect to `topology.get_graph().get_vertices()`
        :return: a reference to this recipe to enable a fluent interface
        """
        self._recipe.separate_vertex(vertex)
        return self

    def change_topology_type(self, new_type):
        """
        Changes the type of the topology to the given type, potentially changing its structural and spatial
        topology reactions.

        :param new_type: the target type
        :return: a reference to this recipe to enable a fluent interface
        """
        self._recipe.change_topology_type(new_type)
        return self

    def _get(self):
        """
        Internal usage only
        :return: returns binding object
        """
        return self._recipe
