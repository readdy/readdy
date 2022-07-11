def get_readdy_plugin_dir():
    import os
    from pathlib import Path
    if "READDY_PLUGIN_DIR" in os.environ.keys():
        return Path(os.environ["READDY_PLUGIN_DIR"]).resolve()
    else:
        return (Path(os.environ["CONDA_PREFIX"]) / 'lib' / 'readdy_plugins').resolve()
