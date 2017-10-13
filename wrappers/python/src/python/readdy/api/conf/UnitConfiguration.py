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
Created on 10.10.17

@author: clonker
"""

from . import ureg as _ureg
from . import Q_ as _Q_


class UnitConfiguration(object):
    def __init__(self, length_unit='nanometer', time_unit='nanosecond', energy_unit='kilojoule/mol',
                 temperature_unit='kelvin'):
        self._length_unit = _ureg.parse_units(length_unit)
        self._time_unit = _ureg.parse_units(time_unit)
        self._energy_unit = _ureg.parse_units(energy_unit)
        self._temperature_unit = _ureg.parse_units(temperature_unit)
        self._boltzmann = _ureg.parse_units('boltzmann_constant')
        self._avogadro = _ureg.parse_units('avogadro_number')

    @property
    def reg(self):
        return _ureg

    @property
    def boltzmann(self):
        return self._boltzmann

    @property
    def avogadro(self):
        return self._avogadro

    @property
    def temperature_unit(self):
        return self._temperature_unit

    @temperature_unit.setter
    def temperature_unit(self, value):
        assert isinstance(value, _ureg.Unit), "temperature unit can only be an instance of unit"
        self._temperature_unit = value

    @property
    def length_unit(self):
        return self._length_unit

    @length_unit.setter
    def length_unit(self, value):
        assert isinstance(value, _ureg.Unit), "length unit can only be an instance of unit"
        self._length_unit = value

    @property
    def time_unit(self):
        return self._time_unit

    @time_unit.setter
    def time_unit(self, value):
        assert isinstance(value, _ureg.Unit), "time unit can only be an instance of unit"
        self._time_unit = value

    @property
    def energy_unit(self):
        return self._energy_unit

    @energy_unit.setter
    def energy_unit(self, value):
        assert isinstance(value, _ureg.Unit), "energy_unit unit can only be an instance of unit"
        self._energy_unit = value

    @property
    def force_constant_unit(self):
        return self.energy_unit / (self.length_unit ** 2)

    @property
    def diffusion_constant_unit(self):
        return (self.length_unit ** 2) / self.time_unit

    @property
    def reaction_rate_unit(self):
        return 1/self.time_unit

    def convert(self, value, target_units):
        if not isinstance(value, _Q_):
            return value
        return value.to(target_units).magnitude


class NoUnitConfiguration(UnitConfiguration):
    def __init__(self):
        super().__init__()
        self._length_unit = 1.
        self._time_unit = 1.
        self._temperature_unit = 1.
        self._energy_unit = 1.

    def convert(self, value, target_units):
        if not isinstance(value, _Q_):
            return value
        return value.magnitude
