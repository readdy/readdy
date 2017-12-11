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


class ReactionRegistry(object):
    def __init__(self, context_reactions, units):
        self._reactions = context_reactions
        self._units = units

    def add_decay(self, name, particle_type, rate):
        """
        Registers a decay reaction

            particle_type -> 0

        Particles of type particle_type will disappear with the given rate

        :param name: label of the reaction
        :param particle_type: the particle type undergoing the reaction
        :param rate: microscopic/intrinsic reaction rate constant [1/time]
        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")
        if not isinstance(particle_type, str):
            raise ValueError("particle_type must be a string")
        rate = self._units.convert(rate, self._units.reaction_rate_unit)
        if not rate > 0.:
            raise ValueError("the reaction rate must be positive")
        self._reactions.add_decay(name, particle_type, rate)

    def add_conversion(self, name, type_from, type_to, rate):
        """
        Registers a conversion reaction

            type_from -> type_to

        Particles of type_from will convert into type_to with the given rate

        :param name: label of the reactions
        :param type_from: particle type of educt
        :param type_to: particle type of product
        :param rate: microscopic/intrinsic reaction rate constant [1/time]
        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")
        if not isinstance(type_from, str):
            raise ValueError("type_from must be a string")
        if not isinstance(type_to, str):
            raise ValueError("type_to must be a string")
        rate = self._units.convert(rate, self._units.reaction_rate_unit)
        if not rate > 0.:
            raise ValueError("the reaction rate must be positive")
        self._reactions.add_conversion(name, type_from, type_to, rate)

    def add_fusion(self, name, type_from1, type_from2, type_to, rate, educt_distance, weight1=0.5, weight2=0.5):
        """
        Registers a fusion/association reaction

            type_from1 +(educt_distance) type_from2 -> type_to

        One particle of type_from1 and another particle of type_from2 fuse to form a particle of type_to.
        This reaction will occur with the given rate if the particles are closer than the educt_distance.
        The product will be placed on the connecting axis of the two particles. The weights determine
        where exactly on this axis the product will be placed. By default this is (0.5, 0.5), which corresponds
        to half the distance of educts. The weights should always add up to 1.

        :param name: label of the reactions
        :param type_from1: particle type of first educt
        :param type_from2: particle type of second educt
        :param type_to: particle type of product
        :param rate: microscopic/intrinsic reaction rate constant [1/time]
        :param educt_distance: maximal distance for the reaction to occur [length]
        :param weight1: value between 0 (product placed at educt1) and 1 (product placed at educt2)
        :param weight2: value between 0 (product placed at educt2) and 1 (product placed at educt1)
        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")
        if not isinstance(type_from1, str):
            raise ValueError("type_from1 must be string")
        if not isinstance(type_from2, str):
            raise ValueError("type_from2 must be string")
        if not isinstance(type_to, str):
            raise ValueError("type_to must be string")
        if not 0. <= weight1 <= 1.:
            raise ValueError("weight1 must be in [0,1]")
        if not 0. <= weight2 <= 1.:
            raise ValueError("weight2 must be in [0,1]")
        rate = self._units.convert(rate, self._units.reaction_rate_unit)
        educt_distance = self._units.convert(educt_distance, self._units.length_unit)
        if not rate > 0.:
            raise ValueError("reaction rate must be positive")
        if not educt_distance > 0.:
            raise ValueError("educt_distance must be positive")
        self._reactions.add_fusion(name, type_from1, type_from2, type_to, rate, educt_distance, weight1, weight2)

    def add_fission(self, name, type_from, type_to1, type_to2, rate, product_distance, weight1=0.5, weight2=0.5):
        """
        Registers a fission/dissociation reaction

            type_from -> type_to1 +(product_distance) type_to2

        A particle of type_from dissociates into two particles, of type_to1 and type_to2 respectively.
        This reaction will occur with the given rate. The products will be placed on an randomly oriented axis.
        Their distance on this axis is at most product_distance. The weights determine where the new particles
        will be with respect to the educt. By default the weights are (0.5, 0.5), which corresponds to each product
        having equal distance from the educt. The weights should always add up to 1.

        :param name: label of the reaction
        :param type_from: particle type of educt
        :param type_to1: particle type of first product
        :param type_to2: particle type of second product
        :param rate: microscopic/intrinsic reaction rate constant [1/time]
        :param product_distance: maximal distance at which products will be placed [length]
        :param weight1: value between 0 (product1 placed at educt) and 1 (product2 placed at educt)
        :param weight2: value between 0 (product2 placed at educt) and 1 (product1 placed at educt)
        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")
        if not isinstance(type_from, str):
            raise ValueError("type_from must be string")
        if not isinstance(type_to1, str):
            raise ValueError("type_to1 must be string")
        if not isinstance(type_to2, str):
            raise ValueError("type_to2 must be string")
        if not 0. <= weight1 <= 1.:
            raise ValueError("weight1 must be in [0,1]")
        if not 0. <= weight2 <= 1.:
            raise ValueError("weight2 must be in [0,1]")
        rate = self._units.convert(rate, self._units.reaction_rate_unit)
        product_distance = self._units.convert(product_distance, self._units.length_unit)
        if not rate > 0.:
            raise ValueError("reaction rate must be positive")
        if not product_distance > 0.:
            raise ValueError("product_distance must be positive")
        self._reactions.add_fission(name, type_from, type_to1, type_to2, rate, product_distance, weight1, weight2)

    def add_enzymatic(self, name, type_catalyst, type_from, type_to, rate, educt_distance):
        """
        Registers an enzymatic/catalytic reaction

            type_from +(educt_distance) type_catalyst -> type_to + type_catalyst

        One particle of type_from can convert into type_to, given that another particle of type_catalyst is within
        educt_distance. Then this reaction will occur with the given rate. The product of type_to will have the
        same position as the educt of type_from.

        :param name: label of the reaction
        :param type_catalyst: particle type of catalyst
        :param type_from: particle type of educt
        :param type_to: particle type of product
        :param rate: microscopic/intrinsic reaction rate constant [1/time]
        :param educt_distance: maximal distance for reaction to occur [length]
        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")
        if not isinstance(type_catalyst, str):
            raise ValueError("type_catalyst must be string")
        if not isinstance(type_from, str):
            raise ValueError("type_from must be string")
        if not isinstance(type_to, str):
            raise ValueError("type_to must be string")
        rate = self._units.convert(rate, self._units.reaction_rate_unit)
        educt_distance = self._units.convert(educt_distance, self._units.length_unit)
        if not rate > 0.:
            raise ValueError("reaction rate must be positive")
        if not educt_distance > 0.:
            raise ValueError("educt_distance must be positive")
        self._reactions.add_enzymatic(name, type_catalyst, type_from, type_to, rate, educt_distance)

    def add(self, descriptor, rate):
        """
        Register a reaction according to the descriptor string. Note that the reaction radius will be interpreted
        in default units.

        Examples:
            A -> 0
            A -> B
            A +(2.0) B -> C [0.5, 0.5]
            C -> A +(2.0) B [0.5, 0.5]
            A +(2.0) C -> B + C

        :param descriptor: the descriptor string, default units
        :param rate: microscopic/intrinsic reaction rate constant [1/time]
        """
        rate = self._units.convert(rate, self._units.reaction_rate_unit)
        self._reactions.add(descriptor, rate)
