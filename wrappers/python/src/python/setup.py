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
