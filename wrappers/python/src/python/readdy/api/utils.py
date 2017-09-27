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
Created on 27.09.17

@author: clonker
"""

import numpy as _np
from readdy._internal.readdybinding.common import Vec as _v3


def vec3_of(value):
    if isinstance(value, _v3):
        return value
    assert isinstance(value, (list, tuple, _np.ndarray, )), "value {} was not list, tuple, or ndarray".format(value)
    if isinstance(value, _np.ndarray):
        if value.squeeze().ndim != 1:
            raise ValueError("Invalid shape for vector!")
        value = value.astype(float).squeeze().tolist()
    value = list(value)
    if len(value) != 3:
        raise ValueError("Length of vector can only be 3 but was {}.".format(len(value)))
    return _v3(*value)
