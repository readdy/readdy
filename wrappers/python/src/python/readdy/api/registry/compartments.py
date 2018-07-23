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
