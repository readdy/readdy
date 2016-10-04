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
