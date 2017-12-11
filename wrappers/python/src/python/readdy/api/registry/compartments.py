# coding=utf-8

# Copyright © 2017 Computational Molecular Biology Group,
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
@author: chrisfroe
"""

from readdy.api.utils import vec3_of as _v3_of

class CompartmentRegistry(object):
    def __init__(self, context_compartments, units):
        self._compartments = context_compartments
        self._units = units

    def add_sphere(self, conversions, name, origin, radius, larger_or_less=False):
        """
        Registers a spherical compartment. The sphere is defined by an origin and a radius.
        Depending on the flag `larger_or_less`, the compartment will be the outside (True) of the sphere
        or the inside (False) of the sphere.

        The characteristic function of the sphere compartment for a position x is

            | x - origin | > radius (`larger_or_less` = True)

        or

            | x - origin | < radius (`larger_or_less` = False)

        If the characteristic function evaluates to True for a particle position, the conversions will be applied.

        :param conversions: dictionary of particle types, converting keys to values
        :param name: label for the compartment
        :param origin: origin of the sphere
        :param radius: radius of the sphere
        :param larger_or_less: determines if the compartment is outside/True or inside/False of the sphere
        """
        origin = self._units.convert(origin, self._units.length_unit)
        radius = self._units.convert(radius, self._units.length_unit)
        if not isinstance(conversions, dict) or len(conversions) == 0:
            raise ValueError("conversions must be a dictionary with at least one entry")
        if not isinstance(name, str):
            raise ValueError("name must be a string")
        if not radius > 0.:
            raise ValueError("radius must be positive")
        if not isinstance(larger_or_less, bool):
            raise ValueError("larger_or_less must be a bool")
        self._compartments.add_sphere(conversions, name, _v3_of(origin), radius, larger_or_less)

    def add_plane(self, conversions, name, normal_coefficients, distance, larger_or_less=True):
        """
        Registers a planar compartment. The plane is defined in Hesse normal form by a normal vector and a distance.
        Depending on the flag `larger_or_less`, the compartment will be where the normal points or on the other side.

        The characteristic function of the plane compartment for a position x is

            x * normal_coefficients - distance > 0 (`larger_or_less` = True)

        or

            x * normal_coefficients - distance < 0 (`larger_or_less` = False)

        If the characteristic function evaluates to True for a particle position, the conversions will be applied.

        :param conversions: dictionary of particle types, converting keys to values
        :param name: label for the compartment
        :param normal_coefficients: normal coefficients of the plane according to Hesse normal form
        :param distance: shorted distance of the plane from the origin (0,0,0) according to Hesse normal form
        :param larger_or_less: determines if the compartment is where the normal points/True or on the other side/False
        """
        normal_coefficients = self._units.convert(normal_coefficients, self._units.length_unit)
        distance = self._units.convert(distance, self._units.length_unit)
        if not isinstance(conversions, dict) or len(conversions) == 0:
            raise ValueError("conversions must be a dictionary with at least one entry")
        if not isinstance(name, str):
            raise ValueError("name must be a string")
        if not distance >= 0.:
            raise ValueError("distance must be non-negative")
        if not isinstance(larger_or_less, bool):
            raise ValueError("larger_or_less must be a bool")
        self._compartments.add_plane(conversions, name, _v3_of(normal_coefficients), distance, larger_or_less)
