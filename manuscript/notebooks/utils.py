import io

import beaker
import dask.diagnostics
import wandb
import fsspec
import numpy as np
import pandas as pd
import xarray as xr

from typing import Sequence, Optional
from upath import UPath


PRECIPITATION = "PRATEsfc"
TIME_MEAN_DIAGNOSTICS_PATH = "time_mean_diagnostics.nc"


def beaker_dataset_id_from_name(name):
    client = beaker.Beaker.from_env()
    experiment_client = beaker.services.ExperimentClient(client)
    return experiment_client.get(name).jobs[-1].result.beaker


def open_beaker_dataset(dataset_id: str, path: str, variables=None) -> xr.Dataset:
    """Given a beaker dataset ID and path within dataset, return an xarray
    dataset.

    Note: dataset must fit in memory. Requires h5netcdf backend.  Adapted from
    ace2-paper repo.
    """
    client = beaker.Beaker.from_env()
    file = client.dataset.get_file(dataset_id, path)
    if variables is None:
        return xr.open_dataset(io.BytesIO(file), engine="h5netcdf").load()
    else:
        return xr.open_dataset(io.BytesIO(file), engine="h5netcdf")[variables].load()


def bounds_to_xesmf_bounds(ds, dim):
    bnds = f"{dim}_bnds"
    bounds = xr.concat(
        [ds[bnds].isel(bnds=0), ds[bnds].isel(bnds=1).isel({dim: -1})], dim=dim
    )
    return bounds.rename({dim: f"{dim}_b"}).drop(f"{dim}_b")


def to_xesmf_grid(ds):
    ds["lon_b"] = bounds_to_xesmf_bounds(ds, "lon")
    ds["lat_b"] = bounds_to_xesmf_bounds(ds, "lat")
    return ds[["lon_b", "lat_b", "lon", "lat"]]


def get_regridder():
    import xesmf as xe

    rename = {
        "grid_xt": "lon",
        "grid_yt": "lat",
        "grid_xt_bnds": "lon_bnds",
        "grid_yt_bnds": "lat_bnds",
    }
    four_degree_grid = xr.open_zarr(
        "gs://vcm-ml-raw-flexible-retention/2024-07-03-C96-SHiELD-SOM/regridded-zarrs/gaussian_grid_45_by_90/1xCO2-ic_0001/fluxes_2d.zarr"
    )
    one_degree_grid = xr.open_zarr(
        "gs://vcm-ml-raw-flexible-retention/2024-07-03-C96-SHiELD-SOM/regridded-zarrs/gaussian_grid_180_by_360/1xCO2-ic_0001/fluxes_2d.zarr"
    )
    four_degree_grid = four_degree_grid.rename(rename)
    one_degree_grid = one_degree_grid.rename(rename)
    one_degree = to_xesmf_grid(one_degree_grid)
    four_degree = to_xesmf_grid(four_degree_grid)
    return xe.Regridder(one_degree, four_degree, method="conservative", periodic=True)


def get_land_fraction(grid: str) -> xr.DataArray:
    # Note the land fraction may not be identical on the coastlines for the C96
    # and C24 data regridded to a 4 degree Gaussian grid.  We have not paid
    # that close attention to non-global averages up to this point though, so
    # for now we will not worry about this detail.
    rename = {
        "grid_xt": "lon",
        "grid_yt": "lat",
        "grid_xt_bnds": "lon_bnds",
        "grid_yt_bnds": "lat_bnds",
    }
    ds = xr.open_zarr(
        f"gs://vcm-ml-raw-flexible-retention/2024-07-03-C96-SHiELD-SOM/regridded-zarrs/{grid}/1xCO2-ic_0001/full_state.zarr"
    )
    ds = ds.rename(rename)
    return ds.isel(time=0).land_fraction.drop_vars("time")


def get_weights(grid: str) -> xr.Dataset:
    land_fraction = get_land_fraction(grid)
    non_land_fraction = 1 - land_fraction
    area_weights = np.cos(np.deg2rad(land_fraction.lat)).broadcast_like(land_fraction)
    return (
        xr.concat(
            [
                area_weights,
                area_weights * land_fraction,
                area_weights * non_land_fraction,
            ],
            dim=pd.Index(["global", "land", "non-land"], name="region"),
        )
        .drop_vars(["lon", "lat"])
        .load()
    )


def rename_dataset(ds):
    generated_variables = [v for v in ds.data_vars if "gen_map-" in v]
    generated = ds[generated_variables]
    rename_generated = {v: v.replace("gen_map-", "") for v in generated_variables}
    generated = generated.rename(rename_generated)
    return generated


def safe_concat(datasets, dim, **kwargs):
    variables = set.intersection(*[set(ds.data_vars) for ds in datasets])
    return xr.concat([ds[variables] for ds in datasets], dim=dim, **kwargs)


def dict_to_dataset(dictionary, dim):
    index = pd.Index(dictionary.keys(), name=dim)
    return xr.concat(dictionary.values(), dim=index)


def open_remote_nc(fs, url):
    data = fs.cat(url)
    f = io.BytesIO(data)
    return xr.open_dataset(f).load()


def scale_precipitation(ds):
    if PRECIPITATION in ds:
        ds[PRECIPITATION] = 86400 * ds[PRECIPITATION]
    return ds


def infer_resolution(ds):
    return round(ds["lon"].diff("lon").isel(lon=0).item())


def regrid_dataset(ds, regridder):
    regridded = regridder(ds)
    return regridded.drop_vars(["lat", "lon"])  # Roundoff error prevents alignment


def load_time_mean_data(
    catalog, models, climates, variables=None, target_resolution=1, regridder=None
):
    cases = catalog[catalog["model"].isin(models) & catalog["forcing"].isin(climates)]
    datasets = {}
    for _, case in cases.iterrows():
        ds = open_beaker_dataset(case["beaker_id"], TIME_MEAN_DIAGNOSTICS_PATH)
        ds = rename_dataset(ds)
        if variables is not None:
            ds = ds[variables]

        resolution = infer_resolution(ds)
        if resolution < target_resolution:
            if regridder is None:
                raise ValueError
            ds = regrid_dataset(ds, regridder)

        key = case["model"], case["forcing"], case["initial_condition"]
        datasets[key] = ds

    index = pd.MultiIndex.from_tuples(
        datasets.keys(), names=("model", "climate", "initial_condition")
    )
    combined = safe_concat(datasets.values(), dim="case")
    result = combined.assign_coords(case=index).unstack("case")
    return scale_precipitation(result)


def load_time_mean_spatial_patterns(
    cases: pd.DataFrame, variables: Optional[Sequence[str]] = None, target_resolution: int = 1, regridder=None
):
    datasets = {}
    for _, case in cases.iterrows():
        ds = open_beaker_dataset(case["beaker_id"], TIME_MEAN_DIAGNOSTICS_PATH)
        ds = rename_dataset(ds)
        if variables is not None:
            ds = ds[variables]

        resolution = infer_resolution(ds)
        if resolution < target_resolution:
            if regridder is None:
                raise ValueError
            ds = regrid_dataset(ds, regridder)

        key = case["model"], case["forcing"], case["initial_condition"]
        datasets[key] = ds

    index = pd.MultiIndex.from_tuples(
        datasets.keys(), names=("model", "climate", "initial_condition")
    )
    combined = safe_concat(datasets.values(), dim="case")
    result = combined.assign_coords(case=index).unstack("case")
    return scale_precipitation(result)


def open_catalog():
    return pd.read_csv("../scripts/inference-run-catalog.csv")


def compute_bias(ds, dim, reference):
    return ds.drop_sel({dim: reference}) - ds.sel({dim: reference})
