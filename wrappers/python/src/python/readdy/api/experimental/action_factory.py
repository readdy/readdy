# coding=utf-8

# Copyright © 2020 Computational Molecular Biology Group,
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

from readdy._internal.readdybinding.api.actions import BreakConfig as _BreakConfig


class BreakConfig(object):
    def __init__(self):
        self._conf = _BreakConfig()

    def add_breakable_pair(self, type1, type2, threshold_energy, rate):
        self._conf.add_breakable_pair(type1, type2, threshold_energy, rate)


# todo some argument checking
class ActionFactory(object):
    def __init__(self, sim):
        self._sim = sim

    def integrator_euler_brownian_dynamics(self, time_step):
        return self._sim.create_action_euler_bd(time_step)

    def calculate_forces(self):
        return self._sim.create_action_calculate_forces()

    def create_neighbor_list(self, interaction_distance):
        return self._sim.create_action_create_neighbor_list(interaction_distance)

    def update_neighbor_list(self):
        return self._sim.create_action_update_neighbor_list()

    def clear_neighbor_list(self):
        return self._sim.create_action_clear_neighbor_list()

    def reaction_handler_uncontrolled_approximation(self, time_step):
        return self._sim.create_action_uncontrolled_approximation(time_step)

    def reaction_handler_gillespie(self, time_step):
        return self._sim.create_action_gillespie(time_step)

    def reaction_handler_detailed_balance(self, time_step):
        return self._sim.create_action_detailed_balance(time_step)

    def topology_reaction_handler(self, time_step):
        return self._sim.create_action_evaluate_topology_reactions(time_step)

    def break_bonds(self, time_step, break_config):
        return self._sim.create_action_break_bonds(time_step, break_config._conf)

    def evaluate_observables(self):
        return self._sim.create_action_evaluate_observables()

    def make_checkpoint(self, base_path, max_n_saves):
        return self._sim.create_action_make_checkpoint(base_path, max_n_saves)
