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
