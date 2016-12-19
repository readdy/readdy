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
Created on 19.12.2016

@author: clonker
"""

import h5py as h5
import numpy as np


def is_readdy_trajectory_file(h5file):
    """
    Method determining if a file is a readdy trajectory file (i.e., contains /readdy/trajectory).
    :param h5file: the filename
    :return: True if is a readdy trajectory file, otherwise False
    """
    with h5.File(h5file) as f:
        return "readdy" in f.keys() and "trajectory" in f["readdy"].keys()


class TrajectoryEntry:
    def __init__(self, type_id, t, particle_id, px, py, pz):
        self._type_id = type_id
        self._t = t
        self._id = particle_id
        self._pos = np.array([px, py, pz], dtype=np.float64)

    @property
    def type_id(self):
        return self._type_id

    @property
    def t(self):
        return self._t

    @property
    def id(self):
        return self._id

    @property
    def position(self):
        return self._pos


def to_trajectory_entries(data):
    if isinstance(data, np.ndarray):
        return [[TrajectoryEntry(x[0], x[1], x[2], x[3], x[4], x[5]) for x in arr] for arr in data]
    else:
        return TrajectoryEntry(data[0], data[1], data[2], data[3], data[4], data[5])


class TrajectoryReader(object):
    @classmethod
    def data_set_path(cls):
        return "/readdy/trajectory"

    def __init__(self, h5file):
        assert is_readdy_trajectory_file(h5file), "the provided file was no readdy trajectory file!"
        with h5.File(h5file) as f:
            ds = f[self.data_set_path()]
            self.n_frames = len(ds)
        self._h5file = h5file

    def __getitem__(self, item):
        with h5.File(self._h5file) as f:
            return to_trajectory_entries(f[self.data_set_path()][item])
