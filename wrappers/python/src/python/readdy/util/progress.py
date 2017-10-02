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
Created on 02.10.17

@author: clonker
"""
from progress_reporter import ProgressReporter as _ProgressReporter


class SimulationProgress(_ProgressReporter):
    def __init__(self, n_callbacks):
        self._progress_register(n_callbacks, description="Simulation progress")

    def callback(self, current_t):
        self._progress_update(1)

    def finish(self):
        self._progress_force_finish()
