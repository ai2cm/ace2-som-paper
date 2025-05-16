import fv3config
import pytest
import yaml

from pathlib import Path
from shield_run.append import (
    compose_simulation_command,
    get_destination_dirs,
    get_fv3config_yml,
    get_n_processes,
    get_segment_number,
    move_output,
)
from shield_run.submit import create


@pytest.fixture
def executable(tmp_path):
    executable_path = tmp_path / "SHiELD.x"
    executable_path.touch()
    return executable_path


@pytest.fixture
def config(tmp_path):
    path = tmp_path / "uncompiled_config.yml"
    uncompiled_config = {
        "base_version": "SHiELD/v0.1",
        "initial_conditions": "/default/initial_conditions",
        "namelist": {
            "fv_core_nml": {"layout": [2, 2]},
            "coupler_nml": {
                "current_date": [2000, 1, 1, 0, 0, 0],
                "days": 0,
                "hours": 1,
                "minutes": 0,
            },
        },
        "diag_table": "/path/to/diag_table",
    }
    with path.open("w") as file:
        yaml.dump(uncompiled_config, file)
    return path


@pytest.fixture()
def destination(config, tmp_path):
    path = tmp_path / "destination"
    create(config, path)
    return path


@pytest.fixture()
def destination_in_flux(destination):
    # Intended for testing file movement after a completed segment
    root = get_destination_dirs(destination)
    restart_staging_dir = root.run_dir / "RESTART"
    restart_staging_dir.mkdir()

    coupler_res = restart_staging_dir / "coupler.res"
    fv_core_res = restart_staging_dir / "fv_core.res.nc"

    diag_table = root.run_dir / "diag_table"
    history_file = root.run_dir / "diagnostics.nc"

    stdout = root.run_dir / "fms.out"
    yml = root.run_dir / "fv3config.yml"

    files = [coupler_res, fv_core_res, diag_table, history_file, stdout, yml]
    for file in files:
        file.touch()

    return destination


def test_create(destination):
    # Check that simulation directories were created
    root = get_destination_dirs(destination)
    assert root.run_dir.is_dir()
    assert root.history_dir.is_dir()
    assert root.restart_dir.is_dir()
    assert root.ascii_dir.is_dir()

    # Check that fv3config.yml exists and is the compiled config
    fv3config_yml = get_fv3config_yml(destination)
    assert fv3config_yml.is_file()
    with fv3config_yml.open("r") as file:
        config = fv3config.load(file)
    assert "gfs_physics_nml" in config["namelist"]


def test_get_n_processes(destination):
    expected = 24
    result = get_n_processes(destination)
    assert result == expected


def test_get_segment_number(destination):
    restart_dir = destination / "restart"
    segment_1 = restart_dir / "2000010101"
    segment_1.mkdir()
    segment_2 = restart_dir / "2000010102"
    segment_2.mkdir()

    expected = 2
    result = get_segment_number(destination)
    assert result == expected


def test_move_output(destination_in_flux):
    root = get_destination_dirs(destination_in_flux)
    start_date = "2000010100"
    end_date = "2000010101"
    move_output(start_date, end_date, root)

    assert (root.restart_dir / end_date / "coupler.res").is_file()
    assert (root.restart_dir / end_date / "fv_core.res.nc").is_file()
    assert (root.history_dir / start_date / "diagnostics.nc").is_file()
    assert (root.ascii_dir / start_date / "fms.out").is_file()
    assert (root.ascii_dir / start_date / "fv3config.yml").is_file()

    assert not (root.history_dir / start_date / "diag_table").is_file()
    assert not (root.ascii_dir / start_date / "diag_table").is_file()


def test_compose_simulation_command(executable):
    n_processes = 24
    expected = [
        "srun",
        "--export=ALL",
        "--ntasks=24",
        "--cpus-per-task=1",
        "./SHiELD.x",
    ]
    result = compose_simulation_command(n_processes, executable)
    assert result == expected
