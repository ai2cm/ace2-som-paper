import dask.diagnostics

from kerchunk_utils import open_kerchunk_json


target_grid = snakemake.wildcards["target_grid"]
chunks_scalar = snakemake.config["target_grid_data"][target_grid]["chunks"]["scalar"]
ds = open_kerchunk_json(snakemake.input.kerchunk_json)
ds = ds.squeeze("scalar_axis").drop_vars("scalar_axis")
ds = ds.chunk(chunks_scalar)

with dask.diagnostics.ProgressBar():
    ds.to_zarr(snakemake.output.store)
