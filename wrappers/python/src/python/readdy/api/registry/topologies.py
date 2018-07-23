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
import readdy._internal.readdybinding.api.top as _top
from readdy._internal.readdybinding.api import AnglePotentialConfiguration as _AnglePotentialConfiguration
from readdy._internal.readdybinding.api import BondedPotentialConfiguration as _BondedPotentialConfiguration
from readdy._internal.readdybinding.api import TorsionPotentialConfiguration as _TorsionPotentialConfiguration
from readdy._internal.readdybinding.api.top import Recipe as _Recipe


class TopologyRegistry(object):
    def __init__(self, context_top_registry, units):
        self._registry = context_top_registry
        self._units = units

    def add_type(self, name):
        """
        Registers a new topology type with a given name. It is used for structural and spatial topology reactions.
        :param name: the name
        """
        self._registry.add_type(name)

    def configure_harmonic_bond(self, type1, type2, force_constant=1., length=1.):
        """
        Configures a harmonic bond between an adjacent pair of particles with types `type1` and `type2`. If the types
        are not equal, the sequence in which they are given does not matter.

        :param type1: type of the first particle
        :param type2: type of the second particle
        :param force_constant: the force constant [energy/length**2]
        :param length: length of the harmonic bond [length]
        """
        force_constant = self._units.convert(force_constant, self._units.force_constant_unit)
        length = self._units.convert(length, self._units.length_unit)
        if force_constant <= 0:
            raise ValueError("The force constant must be strictly positive.")
        if length < 0:
            raise ValueError("The length must be non-negative.")
        cfg = _BondedPotentialConfiguration(force_constant, length, "harmonic")
        self._registry.configure_bond_potential(type1, type2, cfg)

    def configure_harmonic_angle(self, type1, type2, type3, force_constant=1., equilibrium_angle=0.):
        """
        Configures a harmonic angle potential between a particle triple of types `type1` -- `type2` -- `type3`. If the
        first and third type do not coincide, the sequence in which they are given does not matter.

        Angles are internally always expressed in radians.

        :param type1: type of the first particle
        :param type2: type of the second particle
        :param type3: type of the third particle
        :param force_constant: the force constant [energy/angle**2]
        :param equilibrium_angle: the equilibrium angle [angle]
        """
        force_constant = self._units.convert(force_constant, self._units.energy_unit / (self._units.reg.radians ** 2))
        if force_constant <= 0:
            raise ValueError("The force constant must be strictly positive.")
        equilibrium_angle = self._units.convert(equilibrium_angle, self._units.reg.radians)
        cfg = _AnglePotentialConfiguration(force_constant, equilibrium_angle, "harmonic")
        self._registry.configure_angle_potential(type1, type2, type3, cfg)

    def configure_cosine_dihedral(self, type1, type2, type3, type4, force_constant=1., multiplicity=1, phi0=0.):
        """
        Configures a cosine dihedral angle potential between a particle quadruple with types as given in the arguments.
        The sequence in which the types can be reversed, so `(type1, type2, type3, type4)` yields the same potential
        terms as `(type4, type3, type2, type1)`.

        Angles are internally always expressed in radians.

        :param type1: type of the first particle
        :param type2: type of the second particle
        :param type3: type of the third particle
        :param type4: type of the fourth particle
        :param force_constant: the force constant [energy/angle**2]
        :param multiplicity: the multiplicity
        :param phi0: the equilibrium angle [angle]
        """
        force_constant = self._units.convert(force_constant, self._units.energy_unit / (self._units.reg.radians ** 2))
        if force_constant <= 0:
            raise ValueError("The force constant must be strictly positive.")
        phi0 = self._units.convert(phi0, self._units.reg.radians)
        cfg = _TorsionPotentialConfiguration(force_constant, multiplicity, phi0, "cos_dihedral")
        self._registry.configure_torsion_potential(type1, type2, type3, type4, cfg)

    def add_spatial_reaction(self, descriptor, rate=1., radius=1.):
        """
        Adds a spatially dependent topology reaction (i.e., fusion / enzymatic kind). Such a reaction has a rate,
        which determines how likely the reaction occurs per time step as well as a radius, which determines the
        radius which is scanned for potential reaction partners.

        Additionally, a descriptor string has to provided, deciding about the nature of the reaction. It consists out
        of a label and a topology-particle-type reaction equation. There are two different possible types, "enzymatic"
        topology reactions and "fusion" topology reactions. An "enzymatic" spatial topology reaction can be described as

            label: T1(P1) + T2(P2) -> T3(P3) + T4(P4)

        meaning that if a particle of type P1 inside a topology of type T1 is close enough to a particle of type P2
        inside a topology of type T2, the topologies and their respective particles can change their types to T3(P3) and
        T4(P4), potentially enabling different spatial and structural topology reactions as well as locally influencing
        the force field.

        An "fusion" spatial topology reaction can be described as

            label: T1(P1) + T2(P2) -> T3(P3 -- P4)

        indicating that between formerly particles of types P1 and P2 a connection should be formed and their types
        should now be P3 and P4, respectively.
        Additionally, one can specify

            label: T1(P1) + T2(P2) -> T3(P3 -- P4) [self=true]

        meaning that particles within one and the same topology can find themselves as reaction partners.

        :param descriptor: the descriptor string
        :param rate: a fixed rate [1/time]
        :param radius: the radius [length]
        """
        rate = self._units.convert(rate, 1 / self._units.time_unit)
        radius = self._units.convert(radius, self._units.length_unit)
        self._registry.add_spatial_reaction(descriptor, rate, radius)

    def add_structural_reaction(self, topology_type, reaction_function, rate_function,
                                raise_if_invalid=True, expect_connected=False):
        """
        Adds a spatially independent structural topology reaction for a certain topology type. It basically consists
        out of two functions:

        * the reaction function, taking a topology object as input and returning a reaction recipe describing what
          the structural changes are to be applied to the topology
        * the rate function, which takes a topology object as input and returns a corresponding fixed rate

        It should be noted that the rate function will only be evaluated upon changes of the topology, i.e., as rarely
        as possible. The reaction function is evaluated when the actual reaction takes place. Also the rate is expected
        to be returned in terms of the magnitude w.r.t. the default units.

        :param topology_type: the topology type for which this reaction is evaluated
        :param reaction_function: the reaction function, as described above
        :param rate_function: the rate function, as described above
        :param raise_if_invalid: raises an error if the outcome of the reaction function is invalid and set to True,
                                 otherwise it will just roll back to the state of before the reaction and print a
                                 warning into the log
        :param expect_connected: can trigger a raise if set to true and the topology's connectivity graph decayed into
                                 two or more independent components, depending on the value of `raise_if_invalid`.
        """
        fun1, fun2 = _top.ReactionFunction(reaction_function), _top.RateFunction(rate_function)
        reaction = _top.StructuralTopologyReaction(fun1, fun2)
        if raise_if_invalid:
            reaction.raise_if_invalid()
        else:
            reaction.roll_back_if_invalid()
        if expect_connected:
            reaction.expect_connected_after_reaction()
        else:
            reaction.create_child_topologies_after_reaction()
        self._registry.add_structural_reaction(topology_type, reaction)

    def add_topology_dissociation(self, topology_type, bond_breaking_rate):
        """
        Adds a (structural) topology dissociation reaction to a certain topology type, i.e., with a rate of
        `n_edges * bond_breaking_rate` an edge will be removed from the graph, possibly yielding sub-topology instances.
        :param topology_type: the type
        :param bond_breaking_rate: the rate per edge, assuming that every edge corresponds to exactly one bond
        """
        bond_breaking_rate = self._units.convert(bond_breaking_rate, 1 / self._units.time_unit)
        def rate_function(topology):
            return bond_breaking_rate * float(len(topology.get_graph().get_edges()))

        def reaction_function(topology):
            edges = topology.get_graph().get_edges()
            edge = edges[_np.random.randint(0, len(edges))]
            return _Recipe(topology).remove_edge(edge)

        self.add_structural_reaction(topology_type, reaction_function, rate_function, False, False)
