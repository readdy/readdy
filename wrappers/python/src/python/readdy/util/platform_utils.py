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
