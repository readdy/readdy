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

try:
    from itertools import izip
except ImportError:
    izip = zip


def is_readdy_trajectory_file(h5file):
    """
    Method determining if a file is a readdy trajectory file (i.e., contains /readdy/trajectory).
    :param h5file: the filename
    :return: True if is a readdy trajectory file, otherwise False
    """
    with h5.File(h5file) as f:
        return "readdy" in f.keys() and "trajectory" in f["readdy"].keys()


class TrajectoryEntry:
    def __init__(self, type_id, t, particle_id, pos, flavor):
        self._type_id = type_id
        self._t = t
        self._id = particle_id
        self._pos = pos
        self._flavor = flavor

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

    @property
    def flavor(self):
        return self._flavor


def to_trajectory_entries(data, time):
    if isinstance(data, np.ndarray):
        return [[TrajectoryEntry(x["typeId"], t, x["id"], x["pos"], x["flavor"]) for x in d] for (d, t) in
                izip(data, time)]
    else:
        return TrajectoryEntry(data["typeId"], time, data["id"], data["pos"], data["flavor"])


def to_trajectory_entries_flat(data, time, limits, slice):
    times = time[slice]
    entries = []
    for limit, t in izip(limits[slice], times):
        frame = []
        for i in range(limit[0], limit[1]):
            entry = data[i]
            frame.append(TrajectoryEntry(entry["typeId"], t, entry["id"], entry["pos"], entry["flavor"]))
        entries.append(frame)
    return entries


class TrajectoryReader(object):
    @classmethod
    def data_set_path(cls, name):
        if len(name) > 0:
            return "/readdy/trajectory/" + name + "/records"
        else:
            return "/readdy/trajectory/records"

    @classmethod
    def time_information_path(cls, name):
        if len(name) > 0:
            return "/readdy/trajectory/" + name + "/time"
        else:
            return "/readdy/trajectory/time"

    @classmethod
    def group_path(cls, name):
        if len(name) > 0:
            return "/readdy/trajectory/" + name + "/"
        else:
            return "/readdy/trajectory/"

    @classmethod
    def limits_path(cls, name):
        if len(name) > 0:
            return "/readdy/trajectory/" + name + "/limits"
        else:
            return "/readdy/trajectory/limits"

    def __init__(self, h5file, name=""):
        assert is_readdy_trajectory_file(h5file), "the provided file was no readdy trajectory file!"
        self._name = name
        with h5.File(h5file) as f:
            ds = f[self.data_set_path(name)]
            self.n_frames = len(ds)
            self._flat = "limits" in f[self.group_path(name)].keys()
        self._h5file = h5file

    def __getitem__(self, item):
        if not self._flat:
            with h5.File(self._h5file) as f:
                return to_trajectory_entries(f[self.data_set_path(self._name)][item],
                                             f[self.time_information_path(self._name)][item])
        else:
            with h5.File(self._h5file) as f:
                return to_trajectory_entries_flat(f[self.data_set_path(self._name)],
                                                  f[self.time_information_path(self._name)],
                                                  f[self.limits_path(self._name)],
                                                  item)
