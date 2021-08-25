def get_readdy_plugin_dir():
    from pathlib import Path
    return (Path(__file__).parent / '..' / 'readdy_plugins').resolve()
