import fv3config
import pytest

from pathlib import Path
from shield_run.compile_config import compile_config, grid_file_assets, get_resolution


def test_get_resolution():
    config = {"namelist": {"fv_core_nml": {"npx": 49, "npy": 49}}}
    result = get_resolution(config)
    expected = "C48"
    assert result == expected


def test_grid_file_assets():
    config = {"namelist": {"fv_core_nml": {"npx": 49, "npy": 49}}}
    root = "/path"
    result = grid_file_assets(config, root)
    expected = [
        fv3config.get_asset_dict(
            "/path/C48",
            "C48_grid.tile1.nc",
            "INPUT/",
            "C48_grid.tile1.nc",
            copy_method="link",
        ),
        fv3config.get_asset_dict(
            "/path/C48",
            "C48_grid.tile2.nc",
            "INPUT/",
            "C48_grid.tile2.nc",
            copy_method="link",
        ),
        fv3config.get_asset_dict(
            "/path/C48",
            "C48_grid.tile3.nc",
            "INPUT/",
            "C48_grid.tile3.nc",
            copy_method="link",
        ),
        fv3config.get_asset_dict(
            "/path/C48",
            "C48_grid.tile4.nc",
            "INPUT/",
            "C48_grid.tile4.nc",
            copy_method="link",
        ),
        fv3config.get_asset_dict(
            "/path/C48",
            "C48_grid.tile5.nc",
            "INPUT/",
            "C48_grid.tile5.nc",
            copy_method="link",
        ),
        fv3config.get_asset_dict(
            "/path/C48",
            "C48_grid.tile6.nc",
            "INPUT/",
            "C48_grid.tile6.nc",
            copy_method="link",
        ),
        fv3config.get_asset_dict(
            "/path/C48", "C48_mosaic.nc", "INPUT/", "grid_spec.nc", copy_method="link"
        ),
    ]
    assert result == expected


@pytest.mark.parametrize("prescribe_grid", [False, True])
@pytest.mark.parametrize(
    "initial_conditions", [None, Path("/override")], ids=lambda x: str(x)
)
def test_compile_config(prescribe_grid, initial_conditions):
    grid_data_root = Path("/fv3config/root/data/grid_data/v1.2")
    default_initial_conditions = "/default/initial/conditions"
    dummy_asset = fv3config.get_asset_dict("/path", "file.nc", "INPUT/", "file.nc")
    config = {
        "base_version": "SHiELD/v0.1",
        "initial_conditions": default_initial_conditions,
        "namelist": {"fv_core_nml": {"npx": 193, "npy": 193}},
        "patch_files": [dummy_asset],
    }
    if prescribe_grid:
        config["namelist"]["fv_grid_nml"] = {"grid_file": "INPUT/grid_spec.nc"}
    result = compile_config(config, initial_conditions, grid_data_root)

    # Check that initial conditions were overridden if requested
    if initial_conditions is None:
        expected_initial_conditions = default_initial_conditions
    else:
        expected_initial_conditions = str(initial_conditions)
    assert result["initial_conditions"] == expected_initial_conditions

    # Check that config was merged with base config
    assert result["namelist"]["fv_core_nml"]["npx"] == 193  # From minimal config
    assert result["namelist"]["fv_core_nml"]["npy"] == 193  # From minimal config
    assert result["namelist"]["fv_core_nml"]["npz"] == 79  # From base config

    # Check that grid was prescribed or not
    if prescribe_grid:
        grid = grid_file_assets(config, grid_data_root)
        assert result["namelist"]["fv_grid_nml"]["grid_file"] == "INPUT/grid_spec.nc"
        assert result["patch_files"] == [dummy_asset] + grid
    else:
        assert result["patch_files"] == [dummy_asset]
