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
import numpy as np


def scale_const_density(scale_factor):
    factors = ps.get_identity_factors()
    factors["n_particles"] = scale_factor
    factors["box_length"] = np.cbrt(scale_factor)
    return factors


def scale_const_volume(scale_factor):
    factors = ps.get_identity_factors()
    factors["n_particles"] = scale_factor
    factors["box_length"] = 2 * scale_factor


if __name__ == '__main__':
    common.set_logging_level("warn", python_console_out=False)
    n_samples = 5
    n_time_steps = 100
    number_factors = np.logspace(0, 4, n_samples)

    for cll in ["DynamicCLL", "CompactCLL", "ContiguousCLL", "CellDecomposition"]:
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(cll))
        for factor in number_factors:
            print("----------------------------- {} -----------------------------".format(factor))
            factors = scale_const_density(factor)
            s = ps.Collisive("CPU", factors, reaction_scheduler="Gillespie")

            # DynamicCLL, CompactCLL, ContiguousCLL, Adaptive, CellDecomposition
            s.sim.set_kernel_config(json.dumps({"CPU": {
                "neighbor_list": {
                    "cll_radius": 1,
                    "type": "{}".format(cll)
                },
                "thread_config": {
                    "n_threads": -1,
                    "mode": "pool"
                }
            }
            }))
            s.run(n_time_steps)
            print(s.performance())
