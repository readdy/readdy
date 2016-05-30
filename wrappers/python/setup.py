import os
from pprint import pprint

from setuptools import setup, find_packages

from setuptools.command.install import install
from distutils.command.build import build

dyn_lib_extensions = [".so", ".dll", ".lib", ".dylib"]
print_sep = "***" * 15


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

        for curr_dir, curr_subdirs, curr_files in os.walk(os.path.join(file_dir, "readdy2")):
            print("\twalking: %s" % curr_dir)
            for f in curr_files:
                if self.is_dynlib(f):
                    print("\t\tfound dynlib %s" % f)
                    target_files.append(os.path.join(curr_dir, f))
        print("\tdynlibs: %s" % target_files)

        # copy resulting tool to library build folder
        internal_build = os.path.join(self.build_lib, "readdy2", "_internal")
        self.mkpath(internal_build)

        if not self.dry_run:
            for target in target_files:
                self.copy_file(target, internal_build)


class ReaDDyInstall(install):
    def run(self):
        # run original install code
        install.run(self)

        # install libs
        internal_build = os.path.join(self.build_lib, "readdy2", "_internal")
        internal_target = os.path.join(self.install_lib, "readdy2", "_internal")
        print("\t setup.py: installing from %s to %s" % (internal_build, internal_target))
        self.copy_tree(internal_build, internal_target)


metadata = dict(
    name='ReaDDy2',
    version='@READDY_VERSION@',
    package_dir={'': '@CMAKE_CURRENT_SOURCE_DIR@'},
    package_data={'readdy2._internal': ["*"]},
    packages=find_packages(where='@CMAKE_CURRENT_SOURCE_DIR@'),
    cmdclass={'build': ReaDDyBuild, 'install': ReaDDyInstall}
)

if __name__ == '__main__':
    print("%s python setup start %s" % (print_sep, print_sep))
    print("calling setup with metadata:")
    pprint(metadata, indent=2, depth=2, width=2)
    setup(**metadata)
    print("%s python setup end %s" % (print_sep, print_sep))
