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
