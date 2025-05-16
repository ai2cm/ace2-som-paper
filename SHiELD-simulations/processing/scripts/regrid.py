import subprocess

import xarray as xr

from metadata import get_resolution, is_non_metric
from pathlib import Path


def infer_fregrid_input_file(params, wildcards):
    segment_directory = params["segment_directory"]
    tape = wildcards["tape"]
    return segment_directory / tape


def infer_sample_file(params, wildcards):
    input_file = infer_fregrid_input_file(params, wildcards)
    sample_file = f"{input_file}.tile1.nc"
    return sample_file


def get_fregrid_input_mosaic(config, resolution):
    input_mosaic = config["cubed_sphere_grid_data"][resolution]
    return Path(input_mosaic)


def infer_fregrid_input_mosaic(config, params, wildcards):
    sample_file = infer_sample_file(params, wildcards)
    ds = xr.open_dataset(sample_file, chunks={})
    resolution = get_resolution(ds)
    return get_fregrid_input_mosaic(config, resolution)


def get_fregrid_output_mosaic(config, target_grid):
    output_mosaic = config["target_grid_data"][target_grid]["mosaic"]
    return Path(output_mosaic)


def infer_fregrid_output_mosaic(config, wildcards):
    target_grid = wildcards["target_grid"]
    return get_fregrid_output_mosaic(config, target_grid)


def get_fregrid_remap_file(config, source_grid, target_grid):
    remap_filename = f"remap_{source_grid}_to_{target_grid}.nc"
    return Path(config["working_directory"]) / remap_filename


def infer_fregrid_remap_file(config, params, wildcards):
    sample_file = infer_sample_file(params, wildcards)
    ds = xr.open_dataset(sample_file, chunks={})
    source_grid = get_resolution(ds)
    target_grid = wildcards["target_grid"]
    return get_fregrid_remap_file(config, source_grid, target_grid)


def infer_fregrid_scalar_field(params, wildcards):
    sample_file = infer_sample_file(params, wildcards)
    ds = xr.open_dataset(sample_file, chunks={})
    scalar_fields = []
    for var in ds.data_vars:
        if is_non_metric(var):
           scalar_fields.append(var)
    scalar_field = ",".join(scalar_fields)
    return scalar_field


def compose_frenctools_command(command, **kwargs):
    command = [command]
    for key, value in kwargs.items():
        command.append(f"--{key}")
        command.append(str(value))
    return command
        

def fregrid(
    input_mosaic: Path,
    output_mosaic: Path,
    input_file: Path,
    scalar_field: str,
    remap_file: Path,
    output_file: Path,
    interp_method: str = "conserve_order1",
    nthreads: int = 1
):
    command = compose_frenctools_command(
        "fregrid",
        input_mosaic=input_mosaic,
        output_mosaic=output_mosaic,
        input_file=input_file,
        scalar_field=scalar_field,
        interp_method=interp_method,
        remap_file=remap_file,
        nthreads=nthreads,
        output_file=output_file,
    )
    subprocess.run(command, check=True)


if __name__ == "__main__":
    config = snakemake.config
    params = snakemake.params
    wildcards = snakemake.wildcards
    
    input_mosaic = infer_fregrid_input_mosaic(config, params, wildcards)
    output_mosaic = infer_fregrid_output_mosaic(config, wildcards)
    input_file = infer_fregrid_input_file(params, wildcards)
    scalar_field = infer_fregrid_scalar_field(params, wildcards)
    remap_file = infer_fregrid_remap_file(config, params, wildcards)
    output_file = snakemake.output.file

    fregrid(
        input_mosaic,
        output_mosaic,
        input_file,
        scalar_field,
        remap_file,
        output_file
    )
