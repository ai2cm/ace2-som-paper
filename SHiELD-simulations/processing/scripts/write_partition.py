import dask.diagnostics

from xpartition_utils import construct_lazy_dataset

# Note we could consider making chunks_initial resolution dependent, but I
# don't think that's worth the complexity right now.
target_grid = snakemake.wildcards["target_grid"]
chunks_initial = snakemake.config["initial_chunks"]
chunks_2d = snakemake.config["target_grid_data"][target_grid]["chunks"]["2D"]
chunks_3d = snakemake.config["target_grid_data"][target_grid]["chunks"]["3D"]
filter_data = snakemake.config["spherical_harmonic_transform_filter"]["filter_data"]
unfiltered_variables = snakemake.config["spherical_harmonic_transform_filter"].get("unfiltered_variables")

ranks = snakemake.params["ranks"]
rank = int(snakemake.params["rank"])
partition_dims = ["time"]  # TODO: could parameterize this

ds = construct_lazy_dataset(
    snakemake.input.kerchunk_json,
    chunks_initial,
    chunks_2d,
    chunks_3d,
    target_grid,
    filter_data=filter_data,
    unfiltered_variables=unfiltered_variables
)
with dask.diagnostics.ProgressBar():
    ds.partition.write(
        snakemake.input.store,
        ranks,
        partition_dims,
        rank
    )
