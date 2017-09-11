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

class ReactionDiffusionSystem(object):

    def __init__(self, **kwargs):
        self._kbt = kwargs['kbt'] if 'kbt' in kwargs.keys() else 1.0
        self._box_size = kwargs['box_size'] if 'box_size' in kwargs.keys() else 1.0
        self._periodic_boundary = kwargs['periodic_boundary'] if 'periodic_boundary' in kwargs.keys() else 1.0
        self._reactions = []
        self._species = []