def get_readdy_plugin_dir():
    import os
    from pathlib import Path
    if "READDY_PLUGIN_DIR" in os.environ:
        return Path(os.environ["READDY_PLUGIN_DIR"]).resolve()
    elif "CONDA_PREFIX" in os.environ:
        return (Path(os.environ["CONDA_PREFIX"]) / 'lib' / 'readdy_plugins').resolve()
    else:
        return Path(__file__).parent.parent / '_internal' / 'plugins'
