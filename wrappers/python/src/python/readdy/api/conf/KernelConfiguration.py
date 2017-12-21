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


class NOOPKernelConfiguration(object):
    def to_json(self):
        return "{}"


class CPUKernelConfiguration(object):
    def __init__(self):
        self._n_threads = -1
        self._cll_radius = 1

    @property
    def n_threads(self):
        return self._n_threads

    @n_threads.setter
    def n_threads(self, value):
        self._n_threads = value

    @property
    def cell_linked_list_radius(self):
        return self._cll_radius

    @cell_linked_list_radius.setter
    def cell_linked_list_radius(self, value):
        if value <= 0:
            raise ValueError("Only strictly positive cell linked list radii permitted!")
        self._cll_radius = value

    def to_json(self):
        import json
        return json.dumps({"CPU": {
            "neighbor_list": {
                "cll_radius": self.cell_linked_list_radius,
            },
            "thread_config": {
                "n_threads": self.n_threads,
            }
        }
        })


class CPULegacyKernelConfiguration(object):
    def __init__(self):
        self._n_threads = -1
        self._thread_mode = "pool"
        self._neighbor_list = "CompactCLL"
        self._cll_radius = 1

    @property
    def n_threads(self):
        return self._n_threads

    @n_threads.setter
    def n_threads(self, value):
        self._n_threads = value

    @property
    def thread_mode(self):
        return self._thread_mode

    @thread_mode.setter
    def thread_mode(self, value):
        if value in ("pool", "std_thread", "std_async"):
            self._thread_mode = value
        else:
            raise ValueError("Invalid thread mode, select one of \"pool\", \"std_thread\", \"std_async\".")

    @property
    def neighbor_list(self):
        return self._neighbor_list

    @neighbor_list.setter
    def neighbor_list(self, value):
        valid_ones = ["DynamicCLL", "CompactCLL", "ContiguousCLL", "Adaptive", "CellDecomposition"]
        if value in valid_ones:
            self._neighbor_list = value
        else:
            raise ValueError("Invalid neighbor list implementation, select one of {}".format(",".join(valid_ones)))

    @property
    def cell_linked_list_radius(self):
        return self._cll_radius

    @cell_linked_list_radius.setter
    def cell_linked_list_radius(self, value):
        if value <= 0:
            raise ValueError("Only strictly positive cell linked list radii permitted!")
        self._cll_radius = value

    def to_json(self):
        import json
        return json.dumps({"CPU_Legacy": {
            "neighbor_list": {
                "cll_radius": self.cell_linked_list_radius,
                "type": "{}".format(self.neighbor_list)
            },
            "thread_config": {
                "n_threads": self.n_threads,
                "mode": self.thread_mode
            }
        }
        })
