import shutil
import ujson

import xarray as xr

from kerchunk.combine import MultiZarrToZarr
from kerchunk_utils import open_kerchunk_json
from metadata import TIME_DIM


def open_sample_kerchunk_json(single_json_files):
    sample_json_file, *_ = single_json_files
    return open_kerchunk_json(sample_json_file)


def get_identical_dims(single_json_files, concat_dim=TIME_DIM):
    ds = open_sample_kerchunk_json(single_json_files)
    identical_dims = []
    for var in ds.variables:
        if concat_dim not in ds[var].dims:
            identical_dims.append(var)
    return identical_dims


def check_static(single_json_files, concat_dim=TIME_DIM):
    ds = open_sample_kerchunk_json(single_json_files)
    return concat_dim not in ds.dims


single_json_files = snakemake.input
combined_json_file = snakemake.output.json
identical_dims = get_identical_dims(single_json_files)
static = check_static(single_json_files)

if static:
    # Nothing to concatenate; just copy one file as the combined json.
    sample_json_file, *_ = single_json_files
    shutil.copyfile(sample_json_file, combined_json_file)
else:
    mzz = MultiZarrToZarr(
        single_json_files,
        coo_map={TIME_DIM: f"cf:{TIME_DIM}"},
        concat_dims=[TIME_DIM],
        identical_dims=identical_dims
    )
    chunks = mzz.translate()
    output = ujson.dumps(chunks).encode()
    with open(combined_json_file, "wb") as file:
        file.write(output)
