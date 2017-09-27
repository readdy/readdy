# coding=utf-8

# Copyright © 2016 Computational Molecular Biology Group,
#                  Freie Universität Berlin (GER)
#
# This file is part of ReaDDy.
#
# ReaDDy is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General
# Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

"""
Created on 08.09.17

@author: clonker
"""

def _parse_save_args(save_args):
    assert save_args is None or isinstance(save_args, (dict, bool)), \
        "save can only be None, bool, or a dictionary, not {}".format(type(save_args))
    if save_args is None or (isinstance(save_args, bool) and not save_args):
        return None, None
    else:
        assert "chunk_size" in save_args.keys(), "save needs to have a \"chunk_size\" key"
        assert "name" in save_args.keys(), "save needs to have a \"name\" key"
        return save_args["name"], save_args["chunk_size"]

class Observables(object):
    def __init__(self, simulation):
        self._simulation = simulation
        self._sim = self._simulation._simulation
        self._observable_handles = []

    def rdf(self, stride, bin_borders, types_count_from, types_count_to, particle_to_density, callback=None,
            save='default'):
        """
        Calculates and possibly records the radial distribution function.

        :param stride: skip `stride` time steps before evaluating the observable again
        :param bin_borders: the bin borders
        :param types_count_from: type(s) to calculate the distances "from", can be string or list of strings
        :param types_count_to: types to calculate the distances "to", can be string or list of strings
        :param particle_to_density: density of "to" particles, scaling factor in the resulting distribution
        :param callback: callback function that has as argument a tuple with two elements, containing lists of
                         scalars each
        :param save: dictionary containing `name` and `chunk_size` or None to not save the observable to file
        """
        if isinstance(save, str) and save == 'default':
            save = {"name": "rdf", "chunk_size": 100}
        if isinstance(types_count_from, str):
            types_count_from = [types_count_from]
        if isinstance(types_count_to, str):
            types_count_to = [types_count_to]

        assert all([isinstance(x, str) for x in types_count_from]), \
            "types_count_from={} has an invalid type".format(types_count_from)
        assert all([isinstance(x, str) for x in types_count_to]), \
            "types_count_to={} has an invalid type".format(types_count_to)
        handle = self._sim.register_observable_radial_distribution(stride, bin_borders, types_count_from, types_count_to,
                                                                   particle_to_density, callback)
        name, chunk_size = _parse_save_args(save)
        if name is not None and chunk_size is not None:
            self._observable_handles.append((name, chunk_size, handle))

    def reactions(self, stride, callback=None, save='default'):
        """
        Records all reactions.

        :param stride: skip `stride` time steps before evaluating the observable again
        :param callback: callback function that has as argument a list of reaction records
        :param save: dictionary containing `name` and `chunk_size` or None to not save the observable to file
        """
        if isinstance(save, str) and save == 'default':
            save = {"name": "reactions", "chunk_size": 100}
        handle = self._sim.register_observable_reactions(stride, callback)

        name, chunk_size = _parse_save_args(save)
        if name is not None and chunk_size is not None:
            self._observable_handles.append((name, chunk_size, handle))

    def particle_positions(self, stride, types=None, callback=None, save='default'):
        """
        Records particle positions of either all particles or a selection of particle types.

        :param stride: skip `stride` time steps before evaluating the observable again
        :param types: types for which to observe the particle positions, can be None for all types
        :param callback: callback function that has as argument a list of 3d vectors
        :param save: dictionary containing `name` and `chunk_size` or None to not save the observable to file
        """
        if isinstance(save, str) and save == 'default':
            save = {"name": "particle_positions", "chunk_size": 100}

        if types is None:
            types = []
        handle = self._sim.register_observable_particle_positions(stride, types, callback)

        name, chunk_size = _parse_save_args(save)
        if name is not None and chunk_size is not None:
            self._observable_handles.append((name, chunk_size, handle))

    def particles(self, stride, callback=None, save='default'):
        """
        Records all particles in the system.

        :param stride: skip `stride` time steps before evaluating the observable again
        :param callback: callback function that has as argument a list of particle objects
        :param save: dictionary containing `name` and `chunk_size` or None to not save the observable to file
        """
        if isinstance(save, str) and save == 'default':
            save = {"name": "particles", "chunk_size": 100}

        handle = self._sim.register_observable_particles(stride, callback)

        name, chunk_size = _parse_save_args(save)
        if name is not None and chunk_size is not None:
            self._observable_handles.append((name, chunk_size, handle))

    def number_of_particles(self, stride, types=None, callback=None, save='default'):
        """
        Records the number of particles in the system.

        :param stride: skip `stride` time steps before evaluating the observable again
        :param types: types for which to observe the number of particles, can be None for all types
        :param callback: callback function that has as argument a list of particle numbers
        :param save: dictionary containing `name` and `chunk_size` or None to not save the observable to file
        """
        if isinstance(save, str) and save == 'default':
            save = {"name": "n_particles", "chunk_size": 100}

        if types is None:
            types = []
        handle = self._sim.register_observable_n_particles(stride, types, callback)

        name, chunk_size = _parse_save_args(save)
        if name is not None and chunk_size is not None:
            self._observable_handles.append((name, chunk_size, handle))

    def forces(self, stride, types=None, callback=None, save='default'):
        """
        Records the forces acting on particles.

        :param stride: skip `stride` time steps before evaluating the observable again
        :param types: types for which to observe the forces on particles, can be None for all types
        :param callback: callback function that has as argument a list of vectors
        :param save: dictionary containing `name` and `chunk_size` or None to not save the observable to file
        """
        if isinstance(save, str) and save == 'default':
            save = {"name": "forces", "chunk_size": 100}

        if types is None:
            types = []
        handle = self._sim.register_observable_forces(stride, types, callback)

        name, chunk_size = _parse_save_args(save)
        if name is not None and chunk_size is not None:
            self._observable_handles.append((name, chunk_size, handle))

    def reaction_counts(self, stride, callback=None, save='default'):
        """
        Records the occurred number of reactions per registered reaction.

        :param stride: skip `stride` time steps before evaluating the observable again
        :param callback: callback function that has as argument a tuple containing two maps, one for reactions with
                         one educt and one for reactions with two educts
        :param save: dictionary containing `name` and `chunk_size` or None to not save the observable to file
        """
        if isinstance(save, str) and save == 'default':
            save = {"name": "reaction_counts", "chunk_size": 100}

        handle = self._sim.register_observable_reaction_counts(stride, callback)

        name, chunk_size = _parse_save_args(save)
        if name is not None and chunk_size is not None:
            self._observable_handles.append((name, chunk_size, handle))
