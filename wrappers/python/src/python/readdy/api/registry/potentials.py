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
Created on 26.09.17

@author: clonker
"""

import numpy as _np
from readdy._internal.readdybinding.common import Vec as _v3


def _v3_of(value):
    if isinstance(value, _v3):
        return value
    if isinstance(value, _np.ndarray):
        if value.squeeze().ndim != 1:
            raise ValueError("Invalid shape for vector!")
        value = value.astype(float).squeeze().tolist()
    value = list(value)
    if len(value) != 3:
        raise ValueError("Length of vector can only be 3 but was {}.".format(len(value)))
    return _v3(*value)


class PotentialRegistry(object):
    def __init__(self, context_top_registry):
        self._registry = context_top_registry

    def add_box(self, particle_type, force_constant, origin, extent):
        """
        Adds a box potential acting with a harmonic force on particles of the given type once they leave the area
        spanned by the cuboid that has `origin` as its front lower left and `origin+extent` as its back upper right
        vertex, respectively.

        :param particle_type: the particle type for which the potential is registered
        :param force_constant: the force constant
        :param origin: the origin of the box
        :param extent: the extent of the box
        """
        self._registry.add_box(particle_type, force_constant, _v3_of(origin), _v3_of(extent))

    def add_harmonic_repulsion(self, particle_type1, particle_type2, force_constant, interaction_distance):
        """
        Adds harmonic repulsion between particles of type `particle_type1` and `particle_type2`. It is possible to have
        the same type specified for `particle_type1` and `particle_type2`.
        The interaction distance specifies the distance at which particles begin to interact with one another.

        :param particle_type1: first particle type
        :param particle_type2: second particle type
        :param force_constant: the force constant
        :param interaction_distance: the interaction distance
        """
        self._registry.add_harmonic_repulsion(particle_type1, particle_type2, force_constant, interaction_distance)

    def add_weak_interaction_piecewise_harmonic(self, particle_type1, particle_type2, force_constant, desired_distance,
                                                depth, cutoff):
        """
        Adds a weak interaction potential between particles of the specified types. This weak interaction is defined
        by three harmonic potential terms and described by a `desired_distance`, i.e., a distance, at which the
        potential energy is lowest inside the interaction radius, a `depth`, denoting the depth of the potential well,
        and a `cutoff`, denoting the distance at which particles begin to interact.

        :param particle_type1: first particle type
        :param particle_type2: second particle type
        :param force_constant: the force constant
        :param desired_distance: the desired distance, i.e., smallest potential energy
        :param depth: depth of the potential well
        :param cutoff: the cutoff radius
        """
        self._registry.add_weak_interaction_piecewise_harmonic(particle_type1, particle_type2, force_constant,
                                                               desired_distance, depth, cutoff)

    def add_lennard_jones(self, particle_type1, particle_type2, m, n, cutoff, shift, epsilon, sigma):
        """
        Adds a m-n-LJ potential with specified cutoff, epsilon, and sigma. `shift` is bool and denotes whether the
        potential energy should be shifted to bridge the gap at `cutoff`.

        :param particle_type1: first particle type
        :param particle_type2: second particle type
        :param m: first exponent
        :param n: second exponent
        :param cutoff: the cutoff radius
        :param shift: whether to shift the potential energy
        :param epsilon: epsilon value
        :param sigma: sigma value
        """
        self._registry.add_lennard_jones(particle_type1, particle_type2, m, n, cutoff, shift, epsilon, sigma)

    def add_screened_electrostatics(self, particle_type1, particle_type2, electrostatic_strength,
                                    inverse_screening_depth, repulsion_strength, repulsion_distance, exponent, cutoff):
        """
        Adds a screened electrostatics potential between pairs of particles of type `particle_type1` and
        `particle_type2`.

        :param particle_type1:
        :param particle_type2:
        :param electrostatic_strength:
        :param inverse_screening_depth:
        :param repulsion_strength:
        :param repulsion_distance:
        :param exponent:
        :param cutoff:
        """
        self._registry.add_screened_electrostatics(particle_type1, particle_type2, electrostatic_strength,
                                                   inverse_screening_depth, repulsion_strength, repulsion_distance,
                                                   exponent, cutoff)
