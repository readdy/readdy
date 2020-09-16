import os
import platform
import re
import subprocess
import sys
from pathlib import Path

import setuptools

from numpy.distutils.command.build_ext import build_ext
from setuptools._distutils.version import LooseVersion

srcdir = str(Path('.').resolve())


class ReaDDyBindingBuild(build_ext):
    """ from https://github.com/pybind/cmake_example/blob/master/setup.py """

    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        # extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        fullname = self.get_ext_fullname(ext.name)
        ext_filename = os.path.join(self.build_lib,
                                    self.get_ext_filename(fullname))
        # required for auto-detection of auxiliary "native" libs
        # if not ext_filename.endswith(os.path.sep):
        #    ext_filename += os.path.sep

        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + os.path.abspath(os.path.dirname(ext_filename)),
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-DREADDY_BUILD_PYPI_EXT:BOOL=ON',
                      '-DREADDY_BUILD_PYTHON_WRAPPER:BOOL=OFF']

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), ext_filename)]
            if sys.maxsize > 2 ** 32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        print(f"calling {['cmake', srcdir] + cmake_args}")
        subprocess.check_call(['cmake', srcdir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.', '--target', 'readdybinding'] + build_args, cwd=self.build_temp)
        subprocess.check_call(['cmake', '--install', '.', '--component', 'readdybinding'], cwd=self.build_temp)


package_path = Path('.') / 'wrappers' / 'python' / 'src' / 'python'

metadata = dict(name='readdy',
                version='2.0.3',
                author='Moritz Hoffmann',
                author_email='clonker@gmail.com',
                description='An iPRD simulator.',
                long_description='',
                url='https://readdy.github.io',
                install_requires=['numpy', 'cmake', 'pint', 'tqdm', 'h5py'],
                cmdclass=dict(build_ext=ReaDDyBindingBuild),
                zip_safe=False, )


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       )

    config.add_subpackage('readdy', subpackage_path=str((package_path / 'readdy').resolve()))
    config.add_extension('readdy._internal.readdybinding', ['x.c'], runtime_library_dirs="$ORIGIN")
    return config


if __name__ == '__main__':
    print(f'running setup with metadata {metadata}')

    from numpy.distutils.core import setup

    metadata['configuration'] = configuration
    setup(**metadata)
