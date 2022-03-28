def get_readdy_plugin_dir():
    import os
    from pathlib import Path
    if os.environ["READDY_PLUGIN_DIR"]:
        return Path(os.environ["READDY_PLUGIN_DIR"]).resolve()
    else:
        return (Path(os.environ["CONDA_PREFIX"]) / 'lib' / 'readdy_plugins').resolve()
