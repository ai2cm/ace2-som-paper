import argparse
import contextlib
import itertools
import os
import shutil
import subprocess
import sys

import fv3config

from pathlib import Path
from shield_run.dates import get_segment_dates
from typing import Iterable, List, NamedTuple


class SimulationRoot(NamedTuple):
    run_dir: Path
    history_dir: Path
    restart_dir: Path
    ascii_dir: Path


@contextlib.contextmanager
def cwd(path):
    cwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(cwd)


def get_fv3config_yml(destination: Path) -> Path:
    return destination / "fv3config.yml"


def get_n_processes(destination) -> int:
    fv3config_yml = get_fv3config_yml(destination)
    with fv3config_yml.open() as file:
        config = fv3config.load(file)
    return fv3config.config.get_n_processes(config)


def get_destination_dirs(destination) -> SimulationRoot:
    run_dir = destination / "rundir"
    history_dir = destination / "history"
    restart_dir = destination / "restart"
    ascii_dir = destination / "ascii"
    return SimulationRoot(run_dir, history_dir, restart_dir, ascii_dir)


def get_segment_number(destination: Path) -> int:
    restart_dir = destination / "restart"
    return len([p for p in restart_dir.iterdir() if p.is_dir()])


def write_run_directory(executable: Path, config: dict, run_dir: Path) -> None:
    if run_dir.exists():
        shutil.rmtree(run_dir)
        run_dir.mkdir()

    fv3config.write_run_directory(config, run_dir)
    shutil.copy(executable, run_dir)


def move_files(destination: Path, *files: Iterable[Path]) -> None:
    destination.mkdir()
    for file in itertools.chain(*files):
        shutil.move(file, destination)


def move_output(
    start_date: str, end_date: str, simulation_root: SimulationRoot
) -> None:
    run_dir = simulation_root.run_dir
    restart_dir = simulation_root.restart_dir
    history_dir = simulation_root.history_dir
    ascii_dir = simulation_root.ascii_dir

    restart_files = (run_dir / "RESTART").glob("*")
    history_files = run_dir.glob("*.nc")
    ascii_files = run_dir.glob("*.out"), run_dir.glob("*.yml")

    move_files(restart_dir / end_date, restart_files)
    move_files(history_dir / start_date, history_files)
    move_files(ascii_dir / start_date, *ascii_files)


def compose_simulation_command(n_processes: int, executable: Path) -> List[str]:
    executable_file = executable.name
    return [
        "srun",
        "--export=ALL",
        f"--ntasks={n_processes}",
        "--cpus-per-task=1",
        f"./{executable_file}",
    ]


def find(path: str) -> Iterable[Path]:
    return Path(".").rglob("*")


def run_segment(executable: Path, config: dict, run_dir: Path) -> None:
    write_run_directory(executable, config, run_dir)
    with cwd(run_dir):
        manifest = find(".")
        with open("preexisting_files.out", "w") as f:
            for path in manifest:
                print(path, file=f)

        n_processes = fv3config.config.get_n_processes(config)
        command = compose_simulation_command(n_processes, executable)
        with open("fms.out", "w") as file:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            for line in process.stdout:
                for output_stream in [sys.stdout, file]:
                    print(line, end="", file=output_stream)
        return process.wait()


def append(executable: Path, destination: Path) -> None:
    fv3config_yml = get_fv3config_yml(destination)
    with fv3config_yml.open() as file:
        config = fv3config.load(file)

    simulation_root = get_destination_dirs(destination)
    run_dir = simulation_root.run_dir

    segment_number = get_segment_number(destination)
    start_date, end_date = get_segment_dates(config, segment_number)

    if segment_number > 0:
        restart_dir = simulation_root.restart_dir
        initial_conditions = str(restart_dir / start_date)
        config = fv3config.enable_restart(config, initial_conditions)

    exit_code = run_segment(executable, config, run_dir)
    if exit_code != 0:
        raise RuntimeError("Segment crashed; exiting.")

    move_output(start_date, end_date, simulation_root)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("destination", type=Path)
    args, extra_args = parser.parse_known_args()

    append(args.executable, args.destination)
