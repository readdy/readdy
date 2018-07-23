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
Created on 26.09.17

@author: clonker
"""


class NOOPKernelConfiguration(object):
    def to_json(self):
        return "{}"


class CPUKernelConfiguration(object):
    def __init__(self):
        self._n_threads = -1
        self._cll_radius = 1

    @property
    def n_threads(self):
        return self._n_threads

    @n_threads.setter
    def n_threads(self, value):
        self._n_threads = value

    @property
    def cell_linked_list_radius(self):
        return self._cll_radius

    @cell_linked_list_radius.setter
    def cell_linked_list_radius(self, value):
        if value <= 0:
            raise ValueError("Only strictly positive cell linked list radii permitted!")
        self._cll_radius = value

    def to_json(self):
        import json
        return json.dumps({"CPU": {
            "neighbor_list": {
                "cll_radius": self.cell_linked_list_radius,
            },
            "thread_config": {
                "n_threads": self.n_threads,
            }
        }
        })
