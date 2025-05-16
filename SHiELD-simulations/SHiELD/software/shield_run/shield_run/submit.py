import argparse
import contextlib
import os
import subprocess

import fv3config

from pathlib import Path
from shield_run.append import get_destination_dirs, get_fv3config_yml, get_n_processes
from shield_run.compile_config import compile_config, DEFAULT_GRID_DATA_ROOT
from typing import List, Optional


SCRIPT = Path(__file__).parent / "run.sh"


@contextlib.contextmanager
def environment(**environ):
    # https://stackoverflow.com/questions/2059482/temporarily-modify-the-current-processs-environment
    old_environ = dict(os.environ)
    os.environ.update(environ)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(old_environ)


def compose_submission_command(
    executable: Path,
    segments: int,
    time: str,
    destination: Path,
    n_processes: int,
) -> List[str]:
    return [
        "sbatch",
        f"--ntasks={n_processes}",
        "--export=SLURM_NTASKS,SBATCH_TIMELIMIT,CONDA_DEFAULT_ENV",
        f"{SCRIPT}",
        f"{SCRIPT}",
        f"{executable}",
        f"{segments}",
        f"{destination}",
    ]


def create(
    uncompiled_config_path: Path,
    destination: Path,
    initial_conditions: Optional[Path] = None,
    grid_data_root: Path = DEFAULT_GRID_DATA_ROOT,
) -> None:
    # Create directories for storing completed segment output
    destination_dirs = get_destination_dirs(destination)
    for destination_dir in destination_dirs:
        destination_dir.mkdir(parents=True)

    # Open uncompiled config from input path
    with uncompiled_config_path.open() as file:
        uncompiled_config = fv3config.load(file)

    # Compile config and store it in the destination root, similar to the
    # fv3net prognostic run
    fv3config_yml = get_fv3config_yml(destination)
    config = compile_config(uncompiled_config, initial_conditions, grid_data_root)
    with fv3config_yml.open("w") as file:
        fv3config.dump(config, file)


def submit_slurm_script(
    executable: Path,
    segments: int,
    time: str,
    destination: Path,
) -> None:
    n_processes = str(get_n_processes(destination))
    command = compose_submission_command(
        executable,
        segments,
        time,
        destination,
        n_processes,
    )
    with environment(SLURM_NTASKS=n_processes, SBATCH_TIMELIMIT=time):
        subprocess.run(command, check=True)


def submit(
    executable: Path,
    config: Path,
    segments: int,
    time: str,
    destination: Path,
    initial_conditions: Optional[Path] = None,
    grid_data_root: Path = DEFAULT_GRID_DATA_ROOT,
    resume: bool = False,
) -> None:
    if not resume:
        create(config, destination, initial_conditions, grid_data_root)
    submit_slurm_script(
        executable,
        segments,
        time,
        destination,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path, help="Path to the SHiELD executable")
    parser.add_argument("config", type=Path, help="Path to the raw uncompiled config")
    parser.add_argument(
        "segments", type=int, help="Number of segments to run the simulation"
    )
    parser.add_argument("time", type=str, help="Wall clock time per segment")
    parser.add_argument("destination", type=Path, help="Path to the output of the run")
    parser.add_argument(
        "--initial-conditions",
        type=Path,
        default=None,
        help=(
            "Path to the initial conditions for overriding what is in the "
            "config (mainly for initial condition ensembles)."
        ),
    )
    parser.add_argument(
        "--grid-data-root",
        type=Path,
        default=DEFAULT_GRID_DATA_ROOT,
        help=(
            "Path to reference grid data in fv3config structure; used only "
            "appending grid files if requested."
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Whether to resume this simulation after an interruption. "
            "Merely avoids the run creation step."
        ),
    )
    args, extra_args = parser.parse_known_args()

    submit(
        args.executable,
        args.config,
        args.segments,
        args.time,
        args.destination,
        args.initial_conditions,
        args.grid_data_root,
        args.resume,
    )
