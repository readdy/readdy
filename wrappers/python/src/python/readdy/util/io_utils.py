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
Obtain information from a readdy output file, if it contains context info, such as
particle types and reactions. This module shall be standalone and thus not depend on readdy.
It shall only use the specification of the output objects.

@author: chrisfroe
@license: LGPL
"""

import numpy as np
import h5py


def get_particle_types(filename, dset_path="readdy/config/particle_types"):
    """Construct a dictionary from type-strings to type-ids for the particle types
    used in the simulation that created the output file.

    :param filename: the readdy h5 file, containing context info
    :param dset_path: path to dataset within h5 file containing particle type info
    :return: dictionary that maps from type-string to type-id
    """
    result = dict()
    with h5py.File(filename, "r") as f:
        if dset_path in f:
            p_types = f[dset_path]
            for p_type in p_types:
                result[p_type["name"]] = p_type["type_id"]
    return result


def get_particle_types_list(fname, dset_path="readdy/config/particle_types"):
    """Construct a list of type strings for the particle types
    used in the simulation that created the output file.

    :param filename: the readdy h5 file, containing context info
    :param dset_path: path to dataset within h5 file containing particle type info
    :return: list of type strings
    """
    result = []
    with h5py.File(fname, "r") as f:
        if dset_path in f:
            p_types = f[dset_path]
            for p_type in p_types:
                result.append(p_type["name"])
    return result


def get_diffusion_constants(filename, dset_path="readdy/config/particle_types"):
    """Construct a dictionary from type-strings to diffusion constants for the particle types
    used in the simulation that created the output file.

    :param filename: the readdy h5 file, containing context info
    :param dset_path: path to dataset within h5 file containing particle type info
    :return: dictionary that maps from type-string to diffusion constant
    """
    result = dict()
    with h5py.File(filename, "r") as f:
        if dset_path in f:
            p_types = f[dset_path]
            for p_type in p_types:
                result[p_type["name"]] = p_type["diffusion_constant"]
    return result


def get_reactions(filename, dset_path="readdy/config/registered_reactions"):
    """Construct a dictionary where keys are reaction names and value is the corresponding reaction info.

    :param filename: the readdy h5 file, containing context info
    :param dset_path: path to reaction info dataset
    :return: list of reactions
    """
    result = dict()
    with h5py.File(filename, "r") as f:
        if dset_path in f:
            reactions = f[dset_path]
            for r in reactions:
                result[r["name"]] = np.copy(r)
    return result
