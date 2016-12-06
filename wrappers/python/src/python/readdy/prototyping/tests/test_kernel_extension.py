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


import unittest

import readdy._internal.prototyping as pr


class SingleCPUExtension(pr.SingleCPUKernel):
    def __init__(self):
        print("calling super init")
        super(SingleCPUExtension, self).__init__()

    def get_kernel_state_model(self):
        print("calling super get kernel state model")
        super(SingleCPUExtension, self).get_kernel_state_model()

def creator():
    return SingleCPUExtension()


class TestKernelExtension(unittest.TestCase):

    def test_load_prototyping_module(self):
        scpu = SingleCPUExtension()
        program_factory = scpu.get_program_factory()
        addp = program_factory.create_add_particles()
        addp.execute()


if __name__ == '__main__':
    unittest.main()