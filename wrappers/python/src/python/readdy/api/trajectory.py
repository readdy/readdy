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
Created on 28.09.17

@author: clonker
"""

import os as _os
import readdy.util.io_utils as _io_utils
from readdy._internal.readdybinding.common.util import read_trajectory as _read_trajectory


class ReactionInfo:
    def __init__(self, name, index, uuid, n_educts, n_products, rate, educt_distance, product_distance, educt_types,
                 product_types, inverse_types_map):
        self._name = name
        self._index = index
        self._id = uuid
        self._n_educts = n_educts
        self._n_products = n_products
        self._rate = rate
        self._educt_distance = educt_distance
        self._product_distance = product_distance
        self._educt_types = educt_types
        self._product_types = product_types
        self._inverse_types_map = inverse_types_map

    @property
    def name(self):
        """
        Returns the name of the reaction as registered in the reaction diffusion system.
        :return: the name of the reaction
        """
        return self._name

    @property
    def n_educts(self):
        """
        Returns the number of educts in this reaction.
        :return: the number of educts
        """
        return self._n_educts

    @property
    def n_products(self):
        """
        Returns the number of products in this reaction.
        :return: the number of products
        """
        return self._n_products

    @property
    def educt_distance(self):
        """
        Return the educt distance. Only meaningful if n_educts > 1.
        :return: the educt distance
        """
        return self._educt_distance

    @property
    def product_distance(self):
        """
        Returns the product distance. Only meaningful if n_products > 1.
        :return: the product distance
        """
        return self._product_distance

    @property
    def educt_types(self):
        """
        Returns the types that are involved as educts.
        :return: list of length n_educts
        """
        return [self._inverse_types_map[x] for x in self._educt_types[:self.n_educts]]

    @property
    def product_types(self):
        """
        Returns the types that are involved as products.
        :return: list of length n_products
        """
        return [self._inverse_types_map[x] for x in self._product_types[:self.n_products]]

    @property
    def index(self):
        """
        Returns the index of this reaction, can be used for identifying reactions in, e.g., the reactions observable.
        :return: the index
        """
        return self._index

    @property
    def reaction_id(self):
        """
        Returns the global unique id of this reaction.
        :return: the id
        """
        return self._id

    @property
    def rate(self):
        """
        Returns the reaction rate.
        :return: the reaction rate
        """
        return self._rate

    @property
    def type(self):
        """
        Returns the type of this reaction. Can be one of "decay", "conversion", "fission", "fusion", "enzymatic".
        :return: the type of this reaction
        """
        if self.n_educts == 1:
            if self.n_products == 0:
                return "decay"
            elif self.n_products == 1:
                return "conversion"
            elif self.n_products == 2:
                return "fission"
            else:
                raise ValueError("this should not happen, the number of educts was {} and the number "
                                 "of products was {}".format(self.n_educts, self.n_products))
        elif self.n_educts == 2:
            if self.n_products == 1:
                return "fusion"
            elif self.n_products == 2:
                return "enzymatic"
            else:
                raise ValueError("this should not happen, the number of educts was {} and the number "
                                 "of products was {}".format(self.n_educts, self.n_products))
        else:
            raise ValueError("this should not happen, the number of educts was {} and the number "
                             "of products was {}".format(self.n_educts, self.n_products))

    def __str__(self):
        result = ""
        result += self.type + " reaction \"{}\": ".format(self.name)
        result += " + ".join(self.educt_types) + " -> "
        result += " + ".join(self.product_types) + " ["
        result += " rate={} ".format(self.rate)
        if self.n_educts > 1:
            result += " educt_distance={} ".format(self.educt_distance)
        if self.n_products > 1:
            result += " product_distance={} ".format(self.product_distance)
        result += "]"
        return result

    __repr__ = __str__


class Trajectory(object):
    def __init__(self, filename, name=""):
        """
        attempts to open the given trajectory
        :param filename: the file name
        :param name: the trajectory name inside the file, as given in the simulation
        """
        assert _os.path.exists(filename), "The file '{}' did not exist!".format(filename)
        self._filename = filename
        self._name = name
        self._diffusion_constants = _io_utils.get_diffusion_constants(filename)
        self._particle_types = _io_utils.get_particle_types(filename)
        self._reactions_order_1 = []
        inverse_types_map = {v: k for k, v in self.particle_types.items()}
        for reaction in _io_utils.get_reactions_order1(filename):
            info = ReactionInfo(reaction["name"], reaction["index"], reaction["id"], reaction["n_educts"],
                                reaction["n_products"], reaction["rate"], reaction["educt_distance"],
                                reaction["product_distance"], reaction["educt_types"], reaction["product_types"],
                                inverse_types_map)
            self._reactions_order_1.append(info)
            self._reactions_order_2 = []
        for reaction in _io_utils.get_reactions_order2(filename):
            info = ReactionInfo(reaction["name"], reaction["index"], reaction["id"], reaction["n_educts"],
                                reaction["n_products"], reaction["rate"], reaction["educt_distance"],
                                reaction["product_distance"], reaction["educt_types"], reaction["product_types"],
                                inverse_types_map)
            self._reactions_order_2.append(info)

    @property
    def diffusion_constants(self):
        """
        Returns a dictionary of particle type -> diffusion constant
        :return: the diffusion constants
        """
        return self._diffusion_constants

    @property
    def particle_types(self):
        """
        Returns a dictionary of particle type -> particle type ID
        :return: the particle types
        """
        return self._particle_types

    @property
    def reactions_order_1(self):
        """
        Returns a list of `ReactionInfo` objects containing information about each reaction of order 1.
        :return: a list of `ReactionInfo` objects
        """
        return self._reactions_order_1

    @property
    def reactions_order_2(self):
        """
        Returns a list of `ReactionInfo` objects containing information about each reaction of order 2.
        :return: a list of `ReactionInfo` objects
        """
        return self._reactions_order_2

    def read(self):
        """
        Reads the trajectory into memory as a list of lists.

        :return: the trajectory
        """
        return _read_trajectory(self._filename, self._name)
