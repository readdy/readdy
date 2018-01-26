# coding=utf-8

# Copyright © 2018 Computational Molecular Biology Group,
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

import numpy as _np


def calculate_pressure(box_volume, kbt, n_particles, virial):
    if isinstance(n_particles, (list, tuple, _np.ndarray)):
        n_particles = _np.sum(_np.array(n_particles).squeeze())
    VP = n_particles * kbt * _np.identity(3, dtype=virial.dtype) + virial
    return _np.mean(VP.diagonal()) / box_volume
