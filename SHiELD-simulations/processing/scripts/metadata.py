import xarray as xr


TIME_DIM = "time"
X_DIM = "grid_xt"
Y_DIM = "grid_yt"
X_DIM_COARSE = "grid_xt_coarse"
Y_DIM_COARSE = "grid_yt_coarse"
Z_DIM = "pfull"

METRIC_VARIABLES = [
    "area",
    "dx",
    "dy",
    "area_coarse",
    "dx_coarse",
    "dy_coarse",
    "average_T1",
    "average_T2",
    "average_DT",
    "time_bnds",
    "lon",
    "lat",
    "latb",
    "lonb",
    "lon_coarse",
    "lat_coarse",
    "lonb_coarse",
    "latb_coarse",
    "grid_lon",
    "grid_lont",
    "grid_lat",
    "grid_latt",
    "grid_lon_coarse",
    "grid_lont_coarse",
    "grid_lat_coarse",
    "grid_latt_coarse",
    X_DIM,
    Y_DIM,
    X_DIM_COARSE,
    Y_DIM_COARSE,
    TIME_DIM,
    Z_DIM
]
HORIZONTALLY_RESOLVED = (
    [X_DIM, Y_DIM, TIME_DIM],
    [X_DIM_COARSE, Y_DIM_COARSE, TIME_DIM],
)
VERTICALLY_RESOLVED = tuple(
    horizontally_resolved + [Z_DIM]
    for horizontally_resolved in HORIZONTALLY_RESOLVED
)
HORIZONTALLY_RESOLVED_INCLUDE_STATIC = (
    [X_DIM, Y_DIM],
    [X_DIM_COARSE, Y_DIM_COARSE],
)


def is_horizontally_resolved(da: xr.DataArray, include_static=False):
    if include_static:
        reference_dims = HORIZONTALLY_RESOLVED_INCLUDE_STATIC
    else:
        reference_dims = HORIZONTALLY_RESOLVED
    return any(
        all(d in da.dims for d in dims)
        for dims in reference_dims
    )


def is_vertically_resolved(da: xr.DataArray):
    return any(
        all(d in da.dims for d in dims)
        for dims in VERTICALLY_RESOLVED
    )


def is_non_metric(var):
    return var not in METRIC_VARIABLES


def get_resolution(ds):
    if X_DIM in ds.sizes:
        resolution = ds.sizes[X_DIM]
    elif X_DIM_COARSE in ds.sizes:
        resolution = ds.sizes[X_DIM_COARSE]
    else:
        raise KeyError(
            "Could not find a horizontal dimension in dataset to infer the "
            "resolution."
        )
    return f"C{resolution}"


def get_lon_dim(ds):
    for dim in ds.dims:
        if ds[dim].attrs.get("axis") == "X" or ds[dim].attrs.get("cartesian_axis") == "X":
            return dim
    raise ValueError("Could not identify a longitude dimension in the dataset")


def get_lat_dim(ds):
    for dim in ds.dims:
        if ds[dim].attrs.get("axis") == "Y" or ds[dim].attrs.get("cartesian_axis") == "Y":
            return dim
    raise ValueError("Could not identify a latitude dimension in the dataset")
