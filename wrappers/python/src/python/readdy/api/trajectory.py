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
Created on 28.09.17

@author: clonker
"""

import os as _os

import h5py as _h5py
import numpy as _np

from readdy._internal.readdybinding.common.util import read_reaction_observable as _read_reaction_observable
from readdy._internal.readdybinding.common.util import read_trajectory as _read_trajectory
from readdy._internal.readdybinding.common.util import read_topologies_observable as _read_topologies
from readdy.util.observable_utils import calculate_pressure as _calculate_pressure

import readdy.util.io_utils as _io_utils


class ReactionInfo:
    def __init__(self, name, uuid, n_educts, n_products, rate, educt_distance, product_distance, educt_types,
                 product_types, inverse_types_map):
        self._name = str(name)
        self._id = uuid
        self._n_educts = n_educts
        self._n_products = n_products
        self._rate = rate
        self._educt_distance = educt_distance
        self._product_distance = product_distance
        self._educt_types = educt_types
        self._product_types = product_types
        self._inverse_types_map = inverse_types_map

    @property
    def name(self):
        """
        Returns the name of the reaction as registered in the reaction diffusion system.
        :return: the name of the reaction
        """
        return self._name

    @property
    def n_educts(self):
        """
        Returns the number of educts in this reaction.
        :return: the number of educts
        """
        return self._n_educts

    @property
    def n_products(self):
        """
        Returns the number of products in this reaction.
        :return: the number of products
        """
        return self._n_products

    @property
    def educt_distance(self):
        """
        Return the educt distance. Only meaningful if n_educts > 1.
        :return: the educt distance
        """
        return self._educt_distance

    @property
    def product_distance(self):
        """
        Returns the product distance. Only meaningful if n_products > 1.
        :return: the product distance
        """
        return self._product_distance

    @property
    def educt_types(self):
        """
        Returns the types that are involved as educts.
        :return: list of length n_educts
        """
        return [self._inverse_types_map[x] for x in self._educt_types[:self.n_educts]]

    @property
    def product_types(self):
        """
        Returns the types that are involved as products.
        :return: list of length n_products
        """
        return [self._inverse_types_map[x] for x in self._product_types[:self.n_products]]

    @property
    def reaction_id(self):
        """
        Returns the global unique id of this reaction.
        :return: the id
        """
        return self._id

    @property
    def rate(self):
        """
        Returns the reaction rate.
        :return: the reaction rate
        """
        return self._rate

    @property
    def type(self):
        """
        Returns the type of this reaction. Can be one of "decay", "conversion", "fission", "fusion", "enzymatic".
        :return: the type of this reaction
        """
        if self.n_educts == 1:
            if self.n_products == 0:
                return "decay"
            elif self.n_products == 1:
                return "conversion"
            elif self.n_products == 2:
                return "fission"
            else:
                raise ValueError("this should not happen, the number of educts was {} and the number "
                                 "of products was {}".format(self.n_educts, self.n_products))
        elif self.n_educts == 2:
            if self.n_products == 1:
                return "fusion"
            elif self.n_products == 2:
                return "enzymatic"
            else:
                raise ValueError("this should not happen, the number of educts was {} and the number "
                                 "of products was {}".format(self.n_educts, self.n_products))
        else:
            raise ValueError("this should not happen, the number of educts was {} and the number "
                             "of products was {}".format(self.n_educts, self.n_products))

    def __str__(self):
        result = ""
        result += self.type + " reaction \"{}\": ".format(self.name)
        result += " + ".join(self.educt_types) + " -> "
        result += " + ".join(self.product_types) + " ["
        result += " rate={} ".format(self.rate)
        if self.n_educts > 1:
            result += " educt_distance={} ".format(self.educt_distance)
        if self.n_products > 1:
            result += " product_distance={} ".format(self.product_distance)
        result += "]"
        return result

    __repr__ = __str__


class GeneralInformation(object):

    def __init__(self, filename):
        import json

        dsname = "readdy/config/general"
        with _h5py.File(filename, "r") as f:
            if not dsname in f:
                raise ValueError("General information was not recorded in the file!")
            j = json.loads(f[dsname].value)
            self._kbt = j['kbt']
            self._box_volume = j['box_volume']
            self._box_size = _np.array(j['box_size'])
            self._pbc = _np.array(j['pbc'])

    @property
    def kbt(self):
        return self._kbt

    @property
    def box_volume(self):
        return self._box_volume

    @property
    def box_size(self):
        return self._box_size

    @property
    def periodic_boundary_conditions(self):
        return self._pbc


class Trajectory(object):
    def __init__(self, filename, name=""):
        """
        attempts to open the given trajectory
        :param filename: the file name
        :param name: the trajectory name inside the file, as given in the simulation
        """
        assert _os.path.exists(filename), "The file '{}' did not exist!".format(filename)
        self._filename = filename
        self._name = name
        self._diffusion_constants = _io_utils.get_diffusion_constants(filename)
        self._particle_types = _io_utils.get_particle_types(filename)
        self._reactions = []
        self._inverse_types_map = {v: k for k, v in self.particle_types.items()}
        self._general = GeneralInformation(filename)
        for _, reaction in _io_utils.get_reactions(filename).items():
            info = ReactionInfo(reaction["name"], reaction["id"], reaction["n_educts"],
                                reaction["n_products"], reaction["rate"], reaction["educt_distance"],
                                reaction["product_distance"], reaction["educt_types"], reaction["product_types"],
                                self._inverse_types_map)
            self._reactions.append(info)

    def species_name(self, id):
        """
        Retrieves the species' name according to its id as saved in some observables.
        :param id: the id
        :return: the species' name
        """
        return self._inverse_types_map[id]

    @property
    def kbt(self):
        return self._general.kbt

    @property
    def box_volume(self):
        return self._general.box_volume

    @property
    def box_size(self):
        return self._general.box_size

    @property
    def periodic_boundary_conditions(self):
        return self._general.periodic_boundary_conditions

    @property
    def diffusion_constants(self):
        """
        Returns a dictionary of particle type -> diffusion constant
        :return: the diffusion constants
        """
        return self._diffusion_constants

    @property
    def particle_types(self):
        """
        Returns a dictionary of particle type -> particle type ID
        :return: the particle types
        """
        return self._particle_types

    @property
    def reactions(self):
        """
        Returns a list of `ReactionInfo` objects containing information about each reaction.
        :return: a list of `ReactionInfo` objects
        """
        return self._reactions

    def convert_to_xyz(self, xyz_filename=None, generate_tcl=True, tcl_with_grid=False, particle_radii=None):
        """
        Converts this trajectory to a xyz file that can be read into VMD. Assuming the TCL script was generated, the
        trajectory can be visualized by `vmd -e traj.xyz.tcl`.

        :param xyz_filename: the output file name, if None, the filename with an .xyz extension will be used
        :param generate_tcl: generates a tcl script that can be used alongside with the xyz file
        :param tcl_with_grid: enables a grid view inside VMD
        :param particle_radii: map particle radii for visualization purposes, e.g., `{"A": 10., "B": .1}`
        """
        from readdy.api.utils import convert_trajectory_to_xyz as to_xyz
        to_xyz(self._filename, self._name, xyz_filename=xyz_filename, generate_tcl=generate_tcl,
               tcl_with_grid=tcl_with_grid, particle_radii=particle_radii)

    def read(self):
        """
        Reads the trajectory into memory as a list of lists.

        :return: the trajectory
        """
        return _read_trajectory(self._filename, self._name)

    def read_observable_particle_positions(self, data_set_name=""):
        """
        Reads back the output of the particle_positions observable.
        :param data_set_name: The data set name as given in the simulation setup
        :return: a tuple of lists, where the first element contains a list of simulation times and the second element
                 contains a list of (N, 3)-shaped arrays, where N is the number of particles in that time step
        """
        with _h5py.File(self._filename, "r") as f:
            if not "readdy/observables/particle_positions/" + data_set_name in f:
                raise ValueError("The particle positions observable was not recorded in the file or recorded under a "
                                 "different name!")
            group = f["readdy/observables/particle_positions/" + data_set_name]
            time = group["time"][:]
            data = group["data"][:]
            return time, data

    def read_observable_particles(self, data_set_name=""):
        """
        Reads back the output of the particles observable.
        :param data_set_name: The data set name as given in the simulation setup.
        :return: a tuple of lists, where:
                    * the first element contains a list of simulation times
                    * the second element contains a of lists of list of type-ids, which then can be made human-readable
                      by calling `species_name(type_id)`
                    * the third element contains  a list of lists of unique ids for each particle
                    * the fourth element contains a list of lists of particle positions
        """
        with _h5py.File(self._filename, "r") as f:
            group_path = "readdy/observables/particles/" + data_set_name
            if not group_path in f:
                raise ValueError("The particles observable was not recorded in the file or recorded under a different "
                                 "name!")
            group = f[group_path]
            types = group["types"][:]
            ids = group["ids"][:]
            positions = group["positions"][:]
            time = group["time"][:]
            return time, types, ids, positions

    def read_observable_rdf(self, data_set_name="rdf"):
        """
        Reads back the output of the rdf observable.
        :param data_set_name: The data set name as given in the simulation setup.
        :return: a tuple of lists containing (simulation time with shape (T,), bin centers with shape (N, ),
                    distribution value with shape (T, N))
        """
        with _h5py.File(self._filename, "r") as f:
            group_path = "readdy/observables/" + data_set_name
            if not group_path in f:
                raise ValueError("The rdf observable was not recorded in the file or recorded under a different name!")
            group = f[group_path]
            time = group["time"][:]
            bin_centers = group["bin_centers"][:]
            distribution = group["distribution"][:]
            return time, bin_centers, distribution

    def read_observable_number_of_particles(self, data_set_name="n_particles"):
        """
        Reads back the output of the "number of particles" observable.
        :param data_set_name: The data set name as given in the simulation setup.
        :return: a tuple of lists containing the simulation time and a list of lists containing the counts for
                    each specified type
        """
        group_path = "readdy/observables/" + data_set_name
        with _h5py.File(self._filename, "r") as f:
            if not group_path in f:
                raise ValueError("The number of particles observable was not recorded in the file or recorded under a "
                                 "different name!")
            time = f[group_path]["time"][:]
            counts = f[group_path]["data"][:]
            return time, counts

    def read_observable_reactions(self, data_set_name="reactions"):
        """
        Reads back the output of the "reactions" observable
        :param data_set_name: The data set name as given in the simulation setup
        :return: a tuple which contains an array corresponding to the time as first entry and a list of lists containing
                 reaction record objects as second entry
        """
        time = None
        group_path = "readdy/observables/" + data_set_name
        with _h5py.File(self._filename, "r") as f:
            if not group_path in f:
                raise ValueError("The reactions observable was not recorded in the file or recorded under a "
                                 "different name!")
            time = f[group_path]["time"][:]
        return time, _read_reaction_observable(self._filename, data_set_name)

    def read_observable_reaction_counts(self, data_set_name="reaction_counts"):
        """
        Reads back the output of the "reaction_counts" observable
        :param data_set_name: The data set name as given in the simulation setup
        :return:
        """
        group_path = "readdy/observables/" + data_set_name
        with _h5py.File(self._filename, "r") as f:
            if not group_path in f:
                raise ValueError("The reaction counts observable was not recorded in the file or recorded under a "
                                 "different name!")
            time = f[group_path]["time"][:]
            counts = {}
            counts_group = f[group_path]["counts"]
            for reaction in self.reactions:
                if str(reaction.reaction_id) in counts_group:
                    counts[reaction.name] = counts_group[str(reaction.reaction_id)][:]
            return time, counts

    def read_observable_energy(self, data_set_name="energy"):
        """
        Reads back the output of the "energy" observable.
        :param data_set_name: The data set name as given in the simulation setup
        :return: a tuple which contains an array corresponding to the time as first entry and a list of scalars
                 representing the system's potential energy per time step
        """
        group_path = "readdy/observables/" + data_set_name
        with _h5py.File(self._filename, "r") as f:
            if not group_path in f:
                raise ValueError("The energy observable was not recorded in the file or recorded under a "
                                 "different name!")
            time = f[group_path]["time"][:]
            energy = f[group_path]["data"][:]
            return time, energy

    def read_observable_forces(self, data_set_name="forces"):
        """
        Reads back the output of the "forces" observable.
        :param data_set_name: The data set name as given in the simulation setup
        :return: a tuple which contains an array corresponding to the time as first entry and a list of arrays
                 containing the particles' forces per time step
        """
        group_path = "readdy/observables/" + data_set_name
        with _h5py.File(self._filename, "r") as f:
            if not group_path in f:
                raise ValueError("The forces observable was not recorded in the file or recorded under a "
                                 "different name!")
            time = f[group_path]["time"][:]
            forces = f[group_path]["data"][:]
            return time, forces

    def read_observable_topologies(self, data_set_name="topologies"):
        """
        Reads back the output of the "topologies" observable
        :param data_set_name: The data set name as given in the simulation setup
        :return: a tuple which contains an array corresponding to the time as first entry and an array containing
                 lists of topologies per recorded time step
        """
        group_path = "readdy/observables/{}".format(data_set_name)
        with _h5py.File(self._filename, "r") as f:
            if not group_path in f:
                raise ValueError("The topologies observable was not recorded in the file or recorded under a different"
                                 "name!")
        return _read_topologies(self._filename, group_path)

    def read_observable_virial(self, data_set_name="virial"):
        """
        Reads back the output of the "virial" observable.
        :param data_set_name: The data set name as given in the simulation setup
        :return: a tuple which contains an array corresponding to the time as first entry and an array containing the
                 corresponding virial per time step
        """
        group_path = "readdy/observables/{}".format(data_set_name)
        with _h5py.File(self._filename, "r") as f:
            if not group_path in f:
                raise ValueError("The virial observable was not recorded in the file or recorded under a "
                                 "different name!")
            time = f[group_path]["time"][:]
            virial = f[group_path]["data"][:]
            return time, virial

    def read_observable_pressure(self, data_set_name="_pressure"):
        """
        Reads back the output of the "pressure" observable. As the pressure can be computed from the number of particles
        and the virial, this actually reads back the n_particles and virial observables. The data_set_name serves as a
        postfix, where the default value corresponds to the data sets as they are created when using the default
        settings of the observable.

        :param data_set_name: the data set name postfix, default="_pressure",
               yielding "n_particles_pressure" and "virial_pressure", respectively
        :return: a tuple which contains an array corresponding to the time as first entry and an array containing the
                 corresponding pressure per time step
        """
        time_n_particles, n_particles = self.read_observable_number_of_particles("n_particles{}".format(data_set_name))
        time_virial, virial = self.read_observable_virial("virial{}".format(data_set_name))

        if not _np.array_equal(time_n_particles, time_virial):
            raise RuntimeError("For Pressure it is required to know the number of particles and the virial. "
                               "However, these two observables were recorded with different strides.")

        pressure = _np.array([_calculate_pressure(self.box_volume, self.kbt, n_particles[i], virial[i])
                              for i in range(len(time_n_particles))])
        return time_virial, pressure
