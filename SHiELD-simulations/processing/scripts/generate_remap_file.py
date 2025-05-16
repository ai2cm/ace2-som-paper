import subprocess

from pathlib import Path
from regrid import (
    compose_frenctools_command,
    get_fregrid_input_mosaic,
    get_fregrid_output_mosaic,
    get_fregrid_remap_file
)


def generate_remap_file(
    input_mosaic: Path,
    output_mosaic: Path,
    remap_file: Path,
    interp_method: str = "conserve_order1",
    nthreads: int = 1
):
    command = compose_frenctools_command(
        "fregrid",
        input_mosaic=input_mosaic,
        output_mosaic=output_mosaic,
        interp_method=interp_method,
        nthreads=nthreads,
        remap_file=remap_file,
    )
    subprocess.run(command, check=True)


if __name__ == "__main__":
    config = snakemake.config
    wildcards = snakemake.wildcards

    source_grid = wildcards["remap_source_grid"]
    target_grid = wildcards["remap_target_grid"]
    
    input_mosaic = get_fregrid_input_mosaic(config, source_grid)
    output_mosaic = get_fregrid_output_mosaic(config, target_grid)
    remap_file = snakemake.output.remap_file

    generate_remap_file(input_mosaic, output_mosaic, remap_file)
