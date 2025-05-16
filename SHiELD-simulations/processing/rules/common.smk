import xarray as xr

from pathlib import Path

CONFIG = config["config_name"]
GCS_BUCKET = config["gcs_bucket"]
OUTPUT_DIRECTORY = Path(config["output_directory"])
WORKING_DIRECTORY = Path(config["working_directory"])
ARCHIVE_DIRECTORY = Path(config["archive_directory"])
TILES = range(1, 7)
XPARTITION_RANKS = [f"{rank:04d}" for rank in range(config["xpartition_ranks"])]
RUNS = config["runs"].keys()
SOURCE_GRIDS = config["cubed_sphere_grid_data"].keys()
TARGET_GRIDS = config["target_grid_data"].keys()


SEGMENTS = {}
TAPES = {}
for run, run_directory in config["runs"].items():
    glob_string = Path(run_directory) / "history" / "{segment}" / "{tape}.tile1.nc"
    segments, tapes = glob_wildcards(glob_string, followlinks=True)
    SEGMENTS[run] = sorted(list(set(segments)))
    TAPES[run] = sorted(list(set(tapes)))


SCALAR_TAPES = {}
for run, run_directory in config["runs"].items():
    glob_string = Path(run_directory) / "history" / "{segment}" / "{tape}.nc"
    segments, tapes = glob_wildcards(glob_string, followlinks=True)
    SCALAR_TAPES[run] = sorted([tape for tape in set(tapes) if "tile" not in tape])


RESTART_SEGMENTS = {}
for run, run_directory in config["runs"].items():
    glob_string = Path(run_directory) / "restart" / "{restart_segment}" / "fv_core.res.nc"
    segments, = glob_wildcards(glob_string, followlinks=True)
    RESTART_SEGMENTS[run] = sorted(list(set(segments)))


def get_run_directory(wildcards):
    run = wildcards["run"]
    run_directory = config["runs"][run]
    return Path(run_directory)


def get_segment_directory(wildcards):
    segment = wildcards["segment"]
    run_directory = get_run_directory(wildcards)
    return run_directory / "history" / segment


def get_lustre_history_file_path(wildcards):
    segment_directory = get_segment_directory(wildcards)
    tape = wildcards["tape"]
    tile = wildcards["tile"]
    return segment_directory / f"{tape}.tile{tile}.nc"


def get_lustre_scalar_history_file_path(wildcards):
    segment_directory = get_segment_directory(wildcards)
    scalar_tape = wildcards["scalar_tape"]
    return segment_directory / f"{scalar_tape}.nc"


def get_archive_history_file_directory(wildcards):
    run = wildcards["run"]
    segment = wildcards["segment"]
    return ARCHIVE_DIRECTORY / run / "history" / segment


def get_lustre_restart_directory(wildcards):
    run_directory = get_run_directory(wildcards)
    segment = wildcards["restart_segment"]
    return directory(run_directory / "restart" / segment)


def get_archive_restart_file_directory(wildcards):
    run = wildcards["run"]
    return ARCHIVE_DIRECTORY / run / "restart"


def get_lustre_ascii_directory(wildcards):
    run_directory = get_run_directory(wildcards)
    segment = wildcards["segment"]
    return directory(run_directory / "ascii" / segment)


def get_archive_ascii_file_directory(wildcards):
    run = wildcards["run"]
    return ARCHIVE_DIRECTORY / run / "ascii"


def regrid_tape_input(wildcards):
    segment_directory = get_segment_directory(wildcards)
    tape = wildcards["tape"]
    input_file = segment_directory / tape
    pattern = f"{input_file}.tile{{tile}}.nc"
    paths = []
    for tile in TILES:
        path = pattern.format(tile=tile)
        paths.append(path)
    return paths


def kerchunk_regridded_tape_input(wildcards):
    run = wildcards["run"]
    tape = wildcards["tape"]
    target_grid = wildcards["target_grid"]
    segments = SEGMENTS[run]
    return expand(
        rules.kerchunk_regridded_segment_tape.output.json,
        run=run,
        tape=tape,
        target_grid=target_grid,
        segment=segments
    )


def kerchunk_scalar_segment_tape_input(wildcards):
    run_directory = get_run_directory(wildcards)
    scalar_tape = wildcards["scalar_tape"]
    segment = wildcards["segment"]
    return run_directory / "history" / segment / f"{scalar_tape}.nc"


def kerchunk_scalar_tape_input(wildcards):
    run = wildcards["run"]
    scalar_tape = wildcards["scalar_tape"]
    target_grid = wildcards["target_grid"]
    segments = SEGMENTS[run]
    return expand(
        rules.kerchunk_scalar_segment_tape.output.json,
        run=run,
        scalar_tape=scalar_tape,
        segment=segments,
        target_grid=target_grid
    )
