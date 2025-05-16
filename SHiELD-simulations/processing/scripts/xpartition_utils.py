import xarray as xr
import xpartition

from kerchunk_utils import open_kerchunk_json
from metadata import (
    get_lon_dim,
    get_lat_dim,
    is_horizontally_resolved,
    is_vertically_resolved
)
from xtorch_harmonics import roundtrip_filter


def rechunk_if_present(ds, chunks):
    filtered_chunks = {k: v for k, v in chunks.items() if k in ds.dims}
    return ds.chunk(filtered_chunks)


def rechunk(ds, chunks_2d, chunks_3d):
    variables = []
    for v, da in ds.data_vars.items():
        if is_vertically_resolved(da):
            da_chunked = da.chunk(chunks_3d)
        elif is_horizontally_resolved(da):
            da_chunked = da.chunk(chunks_2d)
        else:
            da_chunked = da.chunk()
        variables.append(da_chunked)
    rechunked = xr.merge(variables)
    rechunked.attrs = {}
    return rechunked


def pop_data_vars(ds, variables):
    to_pop = [v for v in variables if v in ds.data_vars]
    to_keep = [v for v in ds.data_vars if v not in variables]
    popped = ds[to_pop]
    kept = ds[to_keep]
    return kept, popped


def get_horizontally_resolved_data_vars(ds):
    names = []
    for name, da in ds.data_vars.items():
        if is_horizontally_resolved(da, include_static=True):
            names.append(name)
    return names


def construct_lazy_dataset(
    kerchunk_json,
    chunks_initial,
    chunks_2d,
    chunks_3d,
    target_grid,
    filter_data=True,
    unfiltered_variables=None,
):
    combined = open_kerchunk_json(kerchunk_json)
    combined = rechunk_if_present(combined, chunks_initial)
    combined = combined.drop_encoding()
    rechunked = rechunk(combined, chunks_2d, chunks_3d)
    lon_dim = get_lon_dim(rechunked)
    lat_dim = get_lat_dim(rechunked)
    regridded_vars = get_horizontally_resolved_data_vars(rechunked)
    formatted_regridded_vars = ", ".join([repr(v) for v in regridded_vars])
    history = (
        f"The following variables were regridded to the {target_grid!r} "
        f"grid using fregrid: {formatted_regridded_vars}."
    )
    if filter_data:
        if unfiltered_variables is None:
            unfiltered_variables = []
        filtered, unfiltered = pop_data_vars(rechunked, unfiltered_variables)
        filtered = roundtrip_filter(filtered, lon_dim=lon_dim, lat_dim=lat_dim)
        result = xr.merge([filtered, unfiltered])
        filtered_vars = get_horizontally_resolved_data_vars(filtered)
        formatted_filtered_vars = ", ".join([repr(v) for v in filtered_vars])
        history += (
            f" The following variables were then filtered with a "
            f"spherical harmonic transform roundtrip assuming a "
            f"legendre-gauss grid: {formatted_filtered_vars}."
        )
        return result.assign_attrs(history=history)
    else:
        return rechunked.assign_attrs(history=history)
