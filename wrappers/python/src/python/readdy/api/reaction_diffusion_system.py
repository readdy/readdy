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

from readdy.api.registry.compartments import Compartments as _Compartments
from readdy._internal.readdybinding.api import Context as _Context
from readdy.api.simulation import Simulation
import numpy as np


class ReactionDiffusionSystem(object):
    def __init__(self):
        self._context = _Context()

    @property
    def kbt(self):
        """
        Returns the thermal energy of the system.
        :return: the thermal energy
        """
        return self._context.kbt

    @kbt.setter
    def kbt(self, value):
        """
        Sets the thermal energy of the system.
        :param value: the new thermal energy
        """
        self._context.kbt = value

    @property
    def box_size(self):
        """
        Retrieves the system's current box size
        :return: the box size
        """
        return self._context.box_size

    @box_size.setter
    def box_size(self, value):
        """
        Sets the system's box size
        :param value: the box size (list or tuple of length 3, numpy scalar array of squeezed shape (3,))
        """
        if isinstance(value, np.ndarray):
            if value.squeeze().ndim != 1:
                raise ValueError("Invalid shape for box size!")
            value = value.astype(float).squeeze().tolist()

        value = list(value)
        if len(value) != 3:
            raise ValueError("Invalid length for box size list! Length can only be 3 but was {}.".format(len(value)))
        self._context.box_size = value

    @property
    def periodic_boundary_conditions(self):
        """
        The periodic boundary conditions of the system.

        If boundaries are set to non-periodic, one has to take care that the particle cannot escape the simulation
        box, e.g., with a box potential.

        :return: a list of length 3 containing boolean values indicating whether x-, y-, or z-boundaries are periodic
                 or not
        """
        return self._context.pbc

    @periodic_boundary_conditions.setter
    def periodic_boundary_conditions(self, value):
        """
        Sets the periodic boundary conditions of the system.
        :param value: List or Tuple or NDarray of length 3 containing boolean values.
        """
        if isinstance(value, np.ndarray):
            if value.squeeze().ndim != 1:
                raise ValueError("Invalid shape for box size!")
            value = value.astype(bool).squeeze().tolist()

        value = list(value)
        if len(value) != 3:
            raise ValueError("Invalid length for box size list! Length can only be 3 but was {}.".format(len(value)))
        self._context.pbc = value

    def calculate_max_cutoff(self):
        """
        Calculates the maximal cutoff distance w.r.t. pair potentials, binary reactions and spatially
        dependent topology reactions.

        This quantity is used for constructing neighbor lists and can have a great impact on the performance.

        :return: The maximal cutoff distance
        """
        self._context.configure(False)
        return self._context.calculate_max_cutoff()

    @property
    def box_volume(self):
        """
        Calculates the box volume.
        :return: the box volume
        """
        return self._context.box_volume()

    @property
    def compartments(self):
        return _Compartments(self._context.compartments)

    def simulation(self, kernel="SingleCPU"):
        """
        Generates a simulation object for this reaction diffusion system configuration.
        :param kernel: The kernel that is used to perform the simulation
        :return: the simulation object
        """
        return Simulation(self, kernel)
