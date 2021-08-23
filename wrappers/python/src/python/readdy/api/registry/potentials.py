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
"""

import numpy as _np
from readdy.api.utils import vec3_of as _v3_of

class PotentialRegistry(object):
    def __init__(self, context_top_registry, context_type_registry, units):
        self._registry = context_top_registry
        self._types = context_type_registry
        self._units = units

    def add_harmonic_geometry(self, particle_type, force_constant, geometry, inclusion):
        r"""
        Adds harmonic inclusion or exclusion based on geometry (accessible via `system.geometry`).

        :param particle_type: The particle type.
        :param force_constant: The force constant.
        :param geometry: Geometry.
        :param inclusion: Whether inclusion (True) or exclusion (False).
        """
        force_constant = self._units.convert(force_constant, self._units.force_constant_unit)
        self._registry.add_harmonic_geometry(particle_type, force_constant, geometry, inclusion)

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

    def add_sphere(self, particle_type, force_constant, origin, radius, inclusion: bool):
        """
        Adds a spherical potential that keeps particles of a certain type restrained to the inside or outside of the
        specified sphere.

        :param particle_type: the particle type
        :param force_constant: strength of the potential [energy/length**2]
        :param origin: origin of the sphere [length]
        :param radius: radius of the sphere [length]
        :param inclusion: if true, the potential will include particles, otherwise exclude them from the volume
        """
        force_constant = self._units.convert(force_constant, self._units.force_constant_unit)
        origin = self._units.convert(origin, self._units.length_unit)
        radius = self._units.convert(radius, self._units.length_unit)
        assert radius > 0, "radius has to be positive"
        self._registry.add_sphere(particle_type, force_constant, _v3_of(origin), radius, inclusion)

    def add_capsule(self, particle_type, force_constant, center, direction, length, radius):
        force_constant = self._units.convert(force_constant, self._units.force_constant_unit)
        center = self._units.convert(center, self._units.length_unit)
        direction = self._units.convert(direction, self._units.length_unit)
        length = self._units.convert(length, self._units.length_unit)
        radius = self._units.convert(radius, self._units.length_unit)
        self._registry.add_capsule(particle_type, force_constant, _v3_of(center), _v3_of(direction), length, radius)

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
        assert _np.abs(height) > 0, "magnitude of height has to be positive"
        assert width > 0, "width has to be positive"
        self._registry.add_spherical_barrier(particle_type, height, width, _v3_of(origin), radius)

    def add_cylinder(self, particle_type, force_constant, origin, normal, radius, inclusion: bool):
        """
        A potential that keeps particles inside or outside a cylindrical volume. Particles penetrating the boundary are
        harmonically pushed back again.

        :param particle_type: the particle type
        :param force_constant: the strength of the confining force [energy/length**2]
        :param origin: any point on the axis of the cylinder [length]
        :param normal: direction of the axis of the cylinder [length]
        :param radius: radius of the cylinder [length]
        :param inclusion: if true, the potential will include particles, otherwise exclude them from the volume
        """
        force_constant = self._units.convert(force_constant, self._units.force_constant_unit)
        origin = self._units.convert(origin, self._units.length_unit)
        normal = self._units.convert(normal, self._units.length_unit)
        radius = self._units.convert(radius, self._units.length_unit)
        assert force_constant > 0, "force_constant has to be positive"
        assert radius > 0, "radius has to be positive"
        self._registry.add_cylinder(particle_type, force_constant, _v3_of(origin), _v3_of(normal), radius, inclusion)

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
