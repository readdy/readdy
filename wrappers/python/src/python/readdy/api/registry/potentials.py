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
from readdy.api.utils import vec3_of as _v3_of

class PotentialRegistry(object):
    def __init__(self, context_top_registry, context_type_registry, units):
        self._registry = context_top_registry
        self._types = context_type_registry
        self._units = units

    def add_box(self, particle_type, force_constant, origin, extent):
        """
        Adds a box potential acting with a harmonic force on particles of the given type once they leave the area
        spanned by the cuboid that has `origin` as its front lower left and `origin+extent` as its back upper right
        vertex, respectively.

        :param particle_type: the particle type for which the potential is registered
        :param force_constant: the force constant [energy/length**2]
        :param origin: the origin of the box [length]
        :param extent: the extent of the box [length]
        """
        force_constant = self._units.convert(force_constant, self._units.force_constant_unit)
        origin = self._units.convert(origin, self._units.length_unit)
        extent = self._units.convert(extent, self._units.length_unit)
        self._registry.add_box(particle_type, force_constant, _v3_of(origin), _v3_of(extent))

    def add_harmonic_repulsion(self, particle_type1, particle_type2, force_constant, interaction_distance):
        """
        Adds harmonic repulsion between particles of type `particle_type1` and `particle_type2`. It is possible to have
        the same type specified for `particle_type1` and `particle_type2`.
        The interaction distance specifies the distance at which particles begin to interact with one another.

        :param particle_type1: first particle type
        :param particle_type2: second particle type
        :param force_constant: the force constant [energy/length**2]
        :param interaction_distance: the interaction distance [length]
        """
        force_constant = self._units.convert(force_constant, self._units.force_constant_unit)
        interaction_distance = self._units.convert(interaction_distance, self._units.length_unit)
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
        :param force_constant: the force constant [energy/length**2]
        :param desired_distance: the desired distance, i.e., smallest potential energy [length]
        :param depth: depth of the potential well [energy]
        :param cutoff: the cutoff radius [length]
        """
        force_constant = self._units.convert(force_constant, self._units.force_constant_unit)
        desired_distance = self._units.convert(desired_distance, self._units.length_unit)
        depth = self._units.convert(depth, self._units.energy_unit)
        cutoff = self._units.convert(cutoff, self._units.length_unit)
        self._registry.add_weak_interaction_piecewise_harmonic(particle_type1, particle_type2, force_constant,
                                                               desired_distance, depth, cutoff)

    def add_lennard_jones(self, particle_type1, particle_type2, epsilon, sigma, cutoff=None, m=12, n=6, shift=True):
        """
        Adds a m-n-LJ potential with specified cutoff, epsilon, and sigma. `shift` is bool and denotes whether the
        potential energy should be shifted to bridge the gap at `cutoff`.

        :param particle_type1: first particle type
        :param particle_type2: second particle type
        :param epsilon: epsilon value [energy]
        :param sigma: sigma value [length]
        :param cutoff: the cutoff radius [length], default value results in 2.5*sigma
        :param m: first exponent, default=12
        :param n: second exponent, default=6
        :param shift: whether to shift the potential energy, default=True
        """
        assert isinstance(shift, bool), "shift can only be bool"
        if cutoff is None:
            cutoff = 2.5*sigma
        cutoff = self._units.convert(cutoff, self._units.length_unit)
        epsilon = self._units.convert(epsilon, self._units.energy_unit)
        sigma = self._units.convert(sigma, self._units.length_unit)
        self._registry.add_lennard_jones(particle_type1, particle_type2, m, n, cutoff, shift, epsilon, sigma)

    def add_screened_electrostatics(self, particle_type1, particle_type2, electrostatic_strength,
                                    inverse_screening_depth, repulsion_strength, repulsion_distance, exponent, cutoff):
        """
        Adds a screened electrostatics potential between pairs of particles of type `particle_type1` and
        `particle_type2`.

        :param particle_type1: first particle type
        :param particle_type2: second particle type
        :param electrostatic_strength: the electrostatic strength [energy * length]
        :param inverse_screening_depth: the inverse screening depth [1 / length]
        :param repulsion_strength: the repulsion strength [energy]
        :param repulsion_distance: the repulsion distance [length]
        :param exponent: the exponent
        :param cutoff: the cutoff radius [length]
        """
        electrostatic_strength = self._units.convert(electrostatic_strength,
                                                     self._units.energy_unit * self._units.length_unit)
        inverse_screening_depth = self._units.convert(inverse_screening_depth, 1/self._units.length_unit)
        repulsion_strength = self._units.convert(repulsion_strength, self._units.energy_unit)
        repulsion_distance = self._units.convert(repulsion_distance, self._units.length_unit)
        cutoff = self._units.convert(cutoff, self._units.length_unit)
        self._registry.add_screened_electrostatics(particle_type1, particle_type2, electrostatic_strength,
                                                   inverse_screening_depth, repulsion_strength, repulsion_distance,
                                                   exponent, cutoff)

    def add_sphere_out(self, particle_type, force_constant, origin, radius):
        """
        Adds a spherical potential that keeps particles of a certain type excluded from the inside of the
        specified sphere.

        :param particle_type: the particle type
        :param force_constant: strength of the potential [energy/length**2]
        :param origin: origin of the sphere [length]
        :param radius: radius of the sphere [length]
        """
        force_constant = self._units.convert(force_constant, self._units.force_constant_unit)
        origin = self._units.convert(origin, self._units.length_unit)
        radius = self._units.convert(radius, self._units.length_unit)
        assert radius > 0, "radius has to be positive"
        self._registry.add_sphere_out(particle_type, force_constant, _v3_of(origin), radius)

    def add_sphere_in(self, particle_type, force_constant, origin, radius):
        """
        Adds a spherical potential that keeps particles of a certain type restrained to the inside of the
        specified sphere.

        :param particle_type: the particle type
        :param force_constant: strength of the potential [energy/length**2]
        :param origin: origin of the sphere [length]
        :param radius: radius of the sphere [length]
        """
        force_constant = self._units.convert(force_constant, self._units.force_constant_unit)
        origin = self._units.convert(origin, self._units.length_unit)
        radius = self._units.convert(radius, self._units.length_unit)
        assert radius > 0, "radius has to be positive"
        self._registry.add_sphere_in(particle_type, force_constant, _v3_of(origin), radius)

    def add_spherical_barrier(self, particle_type, height, width, origin, radius):
        """
        A potential that forms a concentric barrier at a certain radius around a given origin. It is given a height
        (in terms of energy) and a width. Note that the height can also be negative, then this potential acts as
        a 'sticky' sphere. The potential consists of harmonic snippets, such that the energy landscape is continuous
        and differentiable, the force is only continuous and not differentiable.

        :param particle_type: the particle type
        :param height: the height of the barrier [energy]
        :param width: the width of the barrier [length]
        :param origin: the origin of the sphere [length]
        :param radius: the radius of the sphere [length]
        """
        height = self._units.convert(height, self._units.energy_unit)
        width = self._units.convert(width, self._units.length_unit)
        origin = self._units.convert(origin, self._units.length_unit)
        radius = self._units.convert(radius, self._units.length_unit)
        assert radius > 0, "radius has to be positive"
        assert height > 0, "height has to be positive"
        assert width > 0, "width has to be positive"
        self._registry.add_spherical_barrier(particle_type, height, width, _v3_of(origin), radius)

    def add_custom_external(self, particle_type, clazz, *args):
        """
        This method allows to add a custom potential to a simulation. An example of how to define a custom potential
        can be found in the source tree under examples/custom_potential.
        
        :param particle_type: the particle type
        :param clazz: class (not instance) of potential
        :param args: additional arguments that get forwarded into the constructor of the potential
        """
        type_id = self._types.id_of(particle_type)
        instance = clazz(type_id, *args)
        self._registry.add_external_order1(instance)

    def add_custom_pair(self, particle_type1, particle_type2, clazz, *args):
        """
        This method allows to add a custom potential to a simulation. An example of how to define a custom potential
        can be found in the source tree under examples/custom_potential.

        :param particle_type1: first particle type
        :param particle_type2: second particle type
        :param clazz: class (not instance) of potential
        :param args: additional arguments that get forwarded into the constructor of the potential
        """
        type_id1 = self._types.id_of(particle_type1)
        type_id2 = self._types.id_of(particle_type2)
        instance = clazz(type_id1, type_id2, *args)
        self._registry.add_external_order2(instance)
