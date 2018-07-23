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

import os
from pprint import pprint

from setuptools import setup, find_packages

from setuptools.command.install import install
from distutils.command.build import build

dyn_lib_extensions = [".so", ".dll", ".lib", ".dylib"]
print_sep = "***" * 15

__version__ = '@READDY_VERSION@'


class ReaDDyBuild(build):
    @staticmethod
    def is_dynlib(file_path):
        for ext in dyn_lib_extensions:
            if file_path.endswith(ext):
                return True
        return False

    def run(self):
        # super build
        build.run(self)

        # current file
        file_path = os.path.realpath(__file__)
        print("\trealpath: %s" % file_path)
        file_dir = os.path.dirname(__file__)
        print("\tdirname: %s" % file_dir)

        target_files = []

        for curr_dir, curr_subdirs, curr_files in os.walk(os.path.join(file_dir, "readdy")):
            print("\twalking: %s" % curr_dir)
            for f in curr_files:
                if self.is_dynlib(f):
                    print("\t\tfound dynlib %s" % f)
                    target_files.append(os.path.join(curr_dir, f))
        print("\tdynlibs: %s" % target_files)

        # copy resulting tool to library build folder
        internal_build = os.path.join(self.build_lib, "readdy", "_internal")
        self.mkpath(internal_build)

        if not self.dry_run:
            for target in target_files:
                self.copy_file(target, internal_build)


class ReaDDyInstall(install):
    def run(self):
        # run original install code
        install.run(self)

        # install libs
        internal_build = os.path.join(self.build_lib, "readdy", "_internal")
        internal_target = os.path.join(self.install_lib, "readdy", "_internal")
        print("\t setup.py: installing from %s to %s" % (internal_build, internal_target))
        self.copy_tree(internal_build, internal_target)


def get_package_dir():
    return os.path.join('@CMAKE_CURRENT_SOURCE_DIR@', "src", "python")


metadata = dict(
    name='ReaDDy',
    version=__version__[1:] if __version__.startswith("v") else __version__,
    package_dir={'': get_package_dir()},
    package_data={'readdy._internal': ["*"]},
    packages=find_packages(where=get_package_dir()),
    cmdclass={'build': ReaDDyBuild, 'install': ReaDDyInstall}
)

if __name__ == '__main__':
    print("%s python setup start %s" % (print_sep, print_sep))
    print("calling setup with metadata:")
    pprint(metadata, indent=2, width=2)
    setup(**metadata)
    print("%s python setup end %s" % (print_sep, print_sep))
