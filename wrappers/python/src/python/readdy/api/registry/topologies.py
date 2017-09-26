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

from readdy._internal.readdybinding.api import AnglePotentialConfiguration as _AnglePotentialConfiguration
from readdy._internal.readdybinding.api import BondedPotentialConfiguration as _BondedPotentialConfiguration
from readdy._internal.readdybinding.api import TorsionPotentialConfiguration as _TorsionPotentialConfiguration


class TopologyRegistry(object):
    def __init__(self, context_top_registry):
        self._registry = context_top_registry

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
        :param force_constant: the force constant
        :param length: length of the harmonic bond
        """
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

        :param type1: type of the first particle
        :param type2: type of the second particle
        :param type3: type of the third particle
        :param force_constant: the force constant
        :param equilibrium_angle: the equilibrium angle
        """
        if force_constant <= 0:
            raise ValueError("The force constant must be strictly positive.")
        cfg = _AnglePotentialConfiguration(force_constant, equilibrium_angle, "harmonic")
        self._registry.configure_angle_potential(type1, type2, type3, cfg)

    def configure_cosine_dihedral(self, type1, type2, type3, type4, force_constant=1., multiplicity=1, phi0=0.):
        """
        Configures a cosine dihedral angle potential between a particle quadruple with types as given in the arguments.
        The sequence in which the types can be reversed, so `(type1, type2, type3, type4)` yields the same potential
        terms as `(type4, type3, type2, type1)`.

        :param type1: type of the first particle
        :param type2: type of the second particle
        :param type3: type of the third particle
        :param type4: type of the fourth particle
        :param force_constant: the force constant
        :param multiplicity: the multiplicity
        :param phi0: the equilibrium angle
        """
        if force_constant <= 0:
            raise ValueError("The force constant must be strictly positive.")
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
        :param rate: a fixed rate
        :param radius: the radius
        """
        self._registry.add_spatial_reaction(descriptor, rate, radius)
