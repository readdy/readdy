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
@author: chrisfroe
"""

import os
import unittest
import tempfile
import shutil
import contextlib
import numpy as np

import readdy._internal.readdybinding.api as api
import readdy.util.platform_utils as putils
import readdy._internal.readdybinding.common.io as io
import readdy.util.io_utils as ioutils
from readdy.util.testing_utils import ReaDDyTestCase


def get_item(name, collection):
    return next(x for x in collection if x["name"] == name)


class TestIOUtils(ReaDDyTestCase):
    @classmethod
    def setUpClass(cls):
        cls.kernel_provider = api.KernelProvider.get()
        cls.kernel_provider.load_from_dir(putils.get_readdy_plugin_dir())
        cls.dir = tempfile.mkdtemp("test-config-io")
        cls.fname = os.path.join(cls.dir, "test_io_utils.h5")

        sim = api.Simulation("CPU")
        sim.context.particle_types.add("A", 1.)
        sim.context.particle_types.add("B", 2.)
        sim.context.particle_types.add("C", 3.)
        sim.context.reactions.add_conversion("mylabel", "A", "B", .00001)
        sim.context.reactions.add_conversion("A->B", "A", "B", 1.)
        fusion_rate = 0.4
        educt_distance = 0.2
        sim.context.reactions.add_fusion("B+C->A", "B", "C", "A", fusion_rate, educt_distance, .5, .5)
        with contextlib.closing(io.File.create(cls.fname)) as f:
            loop = sim.create_loop(.1)
            loop.write_config_to_file(f)
            loop.run(1)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_particle_types_info(self):
        p_types = ioutils.get_particle_types(self.fname)
        # assuming that type ids are in accordance to order of registration
        np.testing.assert_equal(p_types["A"], 0)
        np.testing.assert_equal(p_types["B"], 1)
        np.testing.assert_equal(p_types["C"], 2)

        diff_constants = ioutils.get_diffusion_constants(self.fname)
        np.testing.assert_equal(diff_constants["A"], 1.)
        np.testing.assert_equal(diff_constants["B"], 2.)
        np.testing.assert_equal(diff_constants["C"], 3.)

    def test_reaction_info(self):
        p_types = ioutils.get_particle_types(self.fname)
        reactions = ioutils.get_reactions(self.fname)
        self.assertEqual(len(reactions), 3)

        self.assertTrue("mylabel" in reactions)
        mylabel = reactions["mylabel"]
        self.assertEqual(mylabel["n_educts"], 1)
        self.assertEqual(mylabel["n_products"], 1)
        self.assertAlmostEqual(mylabel["rate"], 0.00001)
        self.assertEqual(mylabel["educt_types"][0], p_types["A"])
        self.assertEqual(mylabel["product_types"][0], p_types["B"])

        self.assertTrue("A->B" in reactions)
        atob = reactions["A->B"]
        self.assertEqual(atob["n_educts"], 1)
        self.assertEqual(atob["n_products"], 1)
        self.assertEqual(atob["rate"], 1.)
        self.assertEqual(atob["educt_types"][0], p_types["A"])
        self.assertEqual(atob["product_types"][0], p_types["B"])

        self.assertTrue("B+C->A" in reactions)
        fusion = reactions["B+C->A"]
        self.assertEqual(fusion["n_educts"], 2)
        self.assertEqual(fusion["n_products"], 1)
        self.assertAlmostEqual(fusion["rate"], 0.4)
        self.assertAlmostEqual(fusion["educt_distance"], 0.2)
        correct_educts = (fusion["educt_types"][0] == p_types["B"] and fusion["educt_types"][1] == p_types["C"])
        correct_educts = correct_educts or (fusion["educt_types"][1] == p_types["B"] and fusion["educt_types"][0] == p_types["C"])
        self.assertTrue(correct_educts)
        self.assertEqual(fusion["product_types"][0], p_types["A"])

if __name__ == '__main__':
    unittest.main()
