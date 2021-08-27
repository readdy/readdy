def get_readdy_plugin_dir():
    import os
    from pathlib import Path
    return (Path(os.environ["CONDA_PREFIX"]) / 'lib' / 'readdy_plugins').resolve()
