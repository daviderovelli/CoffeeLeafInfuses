import pandas as pd
import yaml
from pathlib import Path

########## Load config.yaml ##########
def load_config(path: str, filepaths: bool=True, **kwargs):

    #load config.yaml as dict
    with open(path, "r") as handle:
        config = yaml.safe_load(handle)

    # Navigate through the config dictionary if kwargs are provided
    for key in kwargs.values():
        if key not in config:
            raise KeyError(f"Key '{key}' not found in the config.yaml")
        config = config[key]

    # Convert string paths to Path objects if filepaths is True
    if filepaths:
        config = {k: Path(v) if isinstance(v, str) else v for k, v in config.items()}

    return config