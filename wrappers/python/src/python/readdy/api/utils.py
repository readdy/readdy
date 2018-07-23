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
Created on 27.09.17

@author: clonker
"""

import numpy as _np
from readdy._internal.readdybinding.common import Vec as _v3
from readdy._internal.readdybinding.common.util import convert_readdyviewer as _convert_readdyviewer
from readdy._internal.readdybinding.common.util import convert_xyz as _convert_xyz


def load_trajectory_to_npy(filename, name="", begin=0, end=None, stride=None):
    """
    Loads the trajectory into four numpy arrays containing:

    * the number of particles per frame (shape=(n_frames,))
    * the positions of the particles for each frame (shape=(n_frames, max_n_particles_per_frame, 3))
    * the type ids of each particle for each frame (shape=(n_frames, max_n_particles_per_frame))
    * the particle unique ids for each frame (shape=(n_frames, max_n_particles_per_frame))

    :param filename: the hdf5 file name
    :param name: the name of the trajectory as it was given when recording it
    :return: the four numpy arrays
    """
    if end is None:
        if stride is None:
            return _convert_readdyviewer(filename, name, begin=begin)
        else:
            return _convert_readdyviewer(filename, name, stride=stride)
    else:
        if stride is None:
            return _convert_readdyviewer(filename, name, begin=begin, end=end)
        else:
            return _convert_readdyviewer(filename, name, stride=stride, begin=begin, end=end)


def convert_trajectory_to_xyz(filename, name="", xyz_filename=None, generate_tcl=True, tcl_with_grid=False,
                              particle_radii=None):
    """
    Converts a hdf5 trajectory to a xyz file that can be read into VMD. Assuming the TCL script was generated, the
    trajectory can be visualized by `vmd -e traj.xyz.tcl`.

    :param filename: the file name of the hdf5 file
    :param name: the name of the trajectory inside the hdf5 file, as given in the simulation configuration
    :param xyz_filename: the output file name, if None, the filename with an .xyz extension will be used
    :param generate_tcl: generates a tcl script that can be used alongside with the xyz file
    :param tcl_with_grid: enables a grid view inside VMD
    :param particle_radii: map particle radii for visualization purposes, e.g., `{"A": 10., "B": .1}`
    """
    if xyz_filename is None:
        xyz_filename = filename + ".xyz"
    if particle_radii is None:
        particle_radii = {}
    _convert_xyz(filename, name, xyz_filename, generate_tcl, tcl_with_grid, particle_radii)


def vec3_of(value):
    """
    Tries to convert the given value to a readdy vec3 object
    :param value: the value
    :return: the vec3 object if successful
    """
    if isinstance(value, _v3):
        return value
    assert isinstance(value, (list, tuple, _np.ndarray,)), "value {} was not list, tuple, or ndarray".format(value)
    if isinstance(value, _np.ndarray):
        if value.squeeze().ndim != 1:
            raise ValueError("Invalid shape for vector!")
        value = value.astype(float).squeeze().tolist()
    value = list(value)
    if len(value) != 3:
        raise ValueError("Length of vector can only be 3 but was {}.".format(len(value)))
    return _v3(*value)
