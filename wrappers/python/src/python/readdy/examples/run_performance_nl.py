#!/usr/bin/env python
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

__author__ = "clonker"

import readdy._internal.readdybinding.common as common

import json

import readdy.examples.performance_scenarios as ps


class ComparePerformance(object):
    def __init__(self, scenario):
        self._scenario = scenario

    @property
    def scenario(self):
        return self._scenario

    def run(self, n_steps):
        common.set_logging_level("debug", python_console_out=False)

        self.scenario.run(n_steps)
        return self.scenario.performance()


if __name__ == '__main__':
    factors = ps.get_identity_factors()



    s = ps.Collisive("CPU", factors, reaction_scheduler="Gillespie")

    s.sim.set_kernel_config(json.dumps({"CPU": {
        "neighbor_list": {
            "cll_radius": 3,
            "type": "CellDecomposition"
        }
    }
    }))
    perf = ComparePerformance(scenario=s)
    times = perf.run(100)
    print(times)
