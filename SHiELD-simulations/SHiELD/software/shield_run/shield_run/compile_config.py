import argparse
import os

import fv3config
import fv3kube
import yaml

from pathlib import Path
from typing import List, Optional


DEFAULT_GRID_DATA_ROOT = "/gpfs/f5/gfdl_w/scratch/Spencer.Clark/reference/2024-01-26-vcm-fv3config/data/grid_data/v1.1"


def grid_file_assets(config: dict, root) -> List[dict]:
    resolution = get_resolution(config)
    files = []
    root = os.path.join(root, resolution)
    for tile in range(1, 7):
        filename = f"{resolution}_grid.tile{tile}.nc"
        asset = fv3config.get_asset_dict(
            root, filename, "INPUT/", filename, copy_method="link"
        )
        files.append(asset)

    filename = f"{resolution}_mosaic.nc"
    mosaic = fv3config.get_asset_dict(
        root, filename, "INPUT/", "grid_spec.nc", copy_method="link"
    )
    files.append(mosaic)
    return files


def get_resolution(config: dict) -> str:
    return f"C{config['namelist']['fv_core_nml']['npx'] - 1}"


def compile_config(
    config: dict,
    initial_conditions: Optional[Path] = None,
    grid_data_root: Path = DEFAULT_GRID_DATA_ROOT,
) -> dict:
    """Compile a minimal config to a full config for use with SHiELD

    Args:
        config: dict
            Minimal config
        initial_conditions: Optional[Path]
            Path to the initial conditions of the simulation
        grid_data_root: Path
            Location of fv3config reference data for grid files (includes the
            version specifier)

    Returns:
        compiled_config: dict
    """
    base_version = config["base_version"]
    base_config = fv3kube.get_base_fv3config(base_version)
    compiled_config = fv3kube.merge_fv3config_overlays(base_config, config)

    # Only override initial conditions from input if it is provided.
    if initial_conditions is not None:
        compiled_config["initial_conditions"] = str(initial_conditions)

    try:
        # Use try-except instead of get due to multiple possibilities for
        # missing keys in nested dict.
        grid_file = compiled_config["namelist"]["fv_grid_nml"]["grid_file"]
    except KeyError:
        grid_file = None

    if grid_file == "INPUT/grid_spec.nc":
        grid_files = grid_file_assets(compiled_config, grid_data_root)
        if "fv_grid_nml" not in compiled_config["namelist"]:
            compiled_config["namelist"]["fv_grid_nml"] = {}
        compiled_config["namelist"]["fv_grid_nml"]["grid_file"] = "INPUT/grid_spec.nc"
        patch_files = compiled_config.get("patch_files", [])
        compiled_config["patch_files"] = patch_files + grid_files

    return compiled_config
