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
