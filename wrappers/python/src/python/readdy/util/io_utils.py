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
Obtain information from a readdy output file, if it contains context info, such as
particle types and reactions. This module shall be standalone and thus not depend on readdy.
It shall only use the specification of the output objects.

@author: chrisfroe
@license: BSD-3
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
