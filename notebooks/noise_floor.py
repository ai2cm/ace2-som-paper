import logging

import dask.diagnostics
import fsspec
import numpy as np
import pandas as pd
import vcm
import xarray as xr

import utils

from upath import UPath


NOISE_FLOOR_DATA_TEMPLATE = "gs://vcm-ml-intermediate/2024-08-09-vertically-resolved-1deg-c96-shield-som-{climate}-dataset-noise-floor"
VERTICAL_COORDINATE_REFERENCE = "gs://vcm-ml-intermediate/2024-07-09-vertically-resolved-1deg-c96-shield-som-ensemble-fme-dataset/1xCO2-ic_0001.zarr"
VERTICALLY_RESOLVED_VARIABLES = [
    "air_temperature",
    "eastward_wind",
    "northward_wind",
    "specific_total_water",
]
LEVELS = 8
UNSTACKED_VERTICALLY_RESOLVED_VARIABLES = []
for variable in VERTICALLY_RESOLVED_VARIABLES:
    levels = [f"{variable}_{n}" for n in range(LEVELS)]
    UNSTACKED_VERTICALLY_RESOLVED_VARIABLES.extend(levels)


def stack_variable(ds, name, levels=8):
    levels = [f"{name}_{n}" for n in range(levels)]
    return ds[levels].to_array(dim="pfull").drop_vars("pfull").rename(name)


def load_vertical_coordinate():
    ds = xr.open_zarr(VERTICAL_COORDINATE_REFERENCE)
    ak = stack_variable(ds, "ak", levels=LEVELS + 1)
    bk = stack_variable(ds, "bk", levels=LEVELS + 1)
    return ak, bk


def compute_delp(surface_pressure):
    ak, bk = load_vertical_coordinate()
    return (ak + surface_pressure * bk).diff("pfull")


def compute_standard_ace_levels():
    ak, bk = load_vertical_coordinate()
    standard_delp = (ak + 100_000 * bk).diff("pfull")
    standard_pressure = vcm.pressure_at_midpoint_log(
        standard_delp, toa_pressure=ak.isel(pfull=0), dim="pfull"
    )
    return standard_delp, standard_pressure


def load_annual_means(climate):
    path = UPath(NOISE_FLOOR_DATA_TEMPLATE.format(climate=climate)) / "annual_means.nc"
    fs, *_ = fsspec.get_fs_token_paths(path)
    return utils.open_remote_nc(fs, str(path)).load()


def load_ensemble_means(climate):
    ds = load_annual_means(climate)
    with dask.diagnostics.ProgressBar():
        result = ds.mean(["sample", "year"]).load()
    return result.rename({"grid_xt": "lon", "grid_yt": "lat"})


def compute_regridded_noise_floor(regridder, climates):
    rmses = {}
    stddevs = {}
    for window_size in [5, 10]:
        rmses[window_size], stddevs[window_size] = combine_climate_stats(
            window_size, 10, 5, regridder, climates
        )
    rmses = xr.concat(rmses.values(), dim=pd.Index(rmses.keys(), name="window_size"))
    stddevs = xr.concat(
        stddevs.values(), dim=pd.Index(stddevs.keys(), name="window_size")
    )
    return rmses, stddevs


def get_stats(
    annual: xr.Dataset,
    years_per_ensemble: int,
    ensemble_members: int,
    window_size: int,
    area: xr.DataArray,
    amip: bool,
):
    """
    Gets the mean and standard deviation of the pattern bias implied by windows
    of a given number of years, from a dataset of annual means.
    """
    n_windows_per_sample = years_per_ensemble // window_size
    rmses = []
    if not amip:
        mean = annual.mean(["sample", "year"])
        annual_bias = annual - mean
    else:
        annual_bias = annual  # bias calculated on each window separately
    for i_w in range(n_windows_per_sample):
        window = annual_bias.isel(
            year=range(i_w * window_size, (i_w + 1) * window_size)
        )
        if amip:
            window = window - window.mean(["sample", "year"])
        bias_maps = window.mean("year")
        rmse = (bias_maps**2).weighted(area).mean(dim=["grid_yt", "grid_xt"]) ** 0.5
        if amip:
            rmse *= (ensemble_members / (ensemble_members - 1.0)) ** 0.5
        else:
            n_windows = n_windows_per_sample * ensemble_members
            rmse *= (n_windows / (n_windows - 1)) ** 0.5
        rmses.append(rmse)
    rmses_da = xr.concat(rmses, "window")
    return rmses_da


def combine_climate_stats(
    window_size: int,
    years_per_ensemble: int,
    ensemble_members: int,
    regridder,
    climates,
    amip: bool = False,
):
    rmses = []
    for climate in climates:
        annual = load_annual_means(climate)
        annual = annual.rename({"grid_xt": "lon", "grid_yt": "lat"})
        regridded = regridder(annual)
        regridded = regridded.rename({"lon": "grid_xt", "lat": "grid_yt"})
        area = np.cos(np.deg2rad(regridded.grid_yt)).broadcast_like(regridded.grid_xt)
        rmse = get_stats(
            regridded, years_per_ensemble, ensemble_members, window_size, area, amip
        )
        rmses.append(rmse)
        del annual
        del regridded
    rmses = xr.concat(rmses, dim="window")
    stddevs = rmses.std(["window", "sample"])
    return rmses.mean(["window", "sample"]), stddevs


def apply_curvefit_result(fit, x, function):
    result = xr.Dataset()
    for variable in fit.data_vars:
        if variable.endswith("_curvefit_coefficients"):
            c = fit[variable].sel(param="c")
            result[variable.replace("_curvefit_coefficients", "")] = function(x, c)
    return result


def extrapolate(ds, window_size=50):
    def function(x, c):
        return c / np.sqrt(x)

    fit = ds.curvefit("window_size", function)
    return apply_curvefit_result(fit, window_size, function)


def compute_regridded_noise_floor_zonal_mean_profile(regridder, climates):
    rmses = {}
    stddevs = {}
    for window_size in [5, 10]:
        rmses[window_size], stddevs[window_size] = (
            combine_climate_stats_zonal_mean_profile(
                window_size, 10, 5, regridder, climates
            )
        )
    rmses = xr.concat(rmses.values(), dim=pd.Index(rmses.keys(), name="window_size"))
    stddevs = xr.concat(
        stddevs.values(), dim=pd.Index(stddevs.keys(), name="window_size")
    )
    return rmses, stddevs


def get_stats_zonal_mean_profile(
    annual: xr.Dataset,
    years_per_ensemble: int,
    ensemble_members: int,
    window_size: int,
    weights: xr.DataArray,
    amip: bool,
):
    """
    Gets the mean and standard deviation of the pattern bias implied by windows
    of a given number of years, from a dataset of annual means.
    """
    n_windows_per_sample = years_per_ensemble // window_size
    rmses = []
    if not amip:
        mean = annual.mean(["sample", "year"])
        annual_bias = annual - mean
    else:
        annual_bias = annual  # bias calculated on each window separately
    for i_w in range(n_windows_per_sample):
        window = annual_bias.isel(
            year=range(i_w * window_size, (i_w + 1) * window_size)
        )
        if amip:
            window = window - window.mean(["sample", "year"])
        bias_maps = window.mean("year")
        rmse = (bias_maps**2).weighted(weights).mean(dim=set(weights.dims) - {"sample"}) ** 0.5
        if amip:
            rmse *= (ensemble_members / (ensemble_members - 1.0)) ** 0.5
        else:
            n_windows = n_windows_per_sample * ensemble_members
            rmse *= (n_windows / (n_windows - 1)) ** 0.5
        rmses.append(rmse)
    rmses_da = xr.concat(rmses, "window")
    return rmses_da


def stack_vertically_resolved_variables(ds):
    arrays = []
    for variable in VERTICALLY_RESOLVED_VARIABLES:
        da = stack_variable(ds, variable, levels=LEVELS)
        arrays.append(da)
    return xr.merge(arrays)


def combine_climate_stats_zonal_mean_profile(
    window_size: int,
    years_per_ensemble: int,
    ensemble_members: int,
    regridder,
    climates,
    amip: bool = False,
):
    rmses = []
    for climate in climates:
        annual = load_annual_means(climate)[
            UNSTACKED_VERTICALLY_RESOLVED_VARIABLES + ["PRESsfc"]
        ]
        annual = annual.rename({"grid_xt": "lon", "grid_yt": "lat"})
        regridded = regridder(annual)
        regridded = regridded.rename({"lon": "grid_xt", "lat": "grid_yt"})

        ak, _ = load_vertical_coordinate()
        stacked_regridded_vertically_resolved = stack_vertically_resolved_variables(
            regridded
        )
        delp = compute_delp(regridded.PRESsfc)
        standard_delp, standard_pressure = compute_standard_ace_levels()
        regridded_interpolated = vcm.interpolate_to_pressure_levels(
            stacked_regridded_vertically_resolved,
            delp,
            levels=standard_pressure.rename({"pfull": "pressure"}),
            dim="pfull",
            ptop=ak.isel(pfull=0),
        )
        regridded_interpolated_zonal_mean = regridded_interpolated.mean("grid_xt")
        weights = standard_delp.rename({"pfull": "pressure"}).where(regridded_interpolated.air_temperature.notnull()).sum("grid_xt") * np.cos(
            np.deg2rad(regridded.grid_yt)
        )
        rmse = get_stats_zonal_mean_profile(
            regridded_interpolated_zonal_mean,
            years_per_ensemble,
            ensemble_members,
            window_size,
            weights,
            amip,
        )
        rmses.append(rmse)
        del annual
        del regridded
    rmses = xr.concat(rmses, dim="window")
    stddevs = rmses.std(["window", "sample"])
    return rmses.mean(["window", "sample"]), stddevs
