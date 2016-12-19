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


from __future__ import print_function
import os


def is_readdy_installed():
    """
    Determines if readdy is installed via checking if the lib/readdy_plugins directory is present.
    :return: True if installed, otherwise False.
    """
    return os.path.exists(os.path.join(get_environment_root_dir(), "readdy", "readdy_plugins"))


def get_readdy_plugin_dir():
    assert is_readdy_installed(), "readdy needs to be installed"
    return os.path.join(get_environment_root_dir(), "readdy", "readdy_plugins")


def get_environment_root_dir():
    """
    Tries to fetch the (conda) environment root dir.
    If the PREFIX environment variable is set, that is what it returns, otherwise it
    looks for the first "lib" directory above os.__file__.

    Raises a ValueError if there is no lib directory in the path.
    :return: The environment root dir
    """
    prefix = os.environ.get('PREFIX')
    if prefix:
        return prefix
    environment_dir = os.__file__
    split = None
    while split != "lib":
        environment_dir, split = os.path.split(environment_dir)
        if not split:
            raise ValueError("Could not find \"lib\" subdir in \"%s\"" % os.__file__)
    return environment_dir
