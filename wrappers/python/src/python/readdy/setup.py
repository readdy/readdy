# numpy distuil setup file, only used for pypi build

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('readdy', parent_package, top_path)
    config.add_subpackage('_internal')
    config.add_subpackage('api')
    config.add_subpackage('api.conf')
    config.add_subpackage('api.experimental')
    config.add_subpackage('api.registry')
    config.add_subpackage('tests')
    config.add_subpackage('util')
    return config
