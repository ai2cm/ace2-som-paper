import xarray as xr


def open_kerchunk_json(path):
    return xr.open_dataset(
        "reference://",
        engine="zarr",
        backend_kwargs={
            "storage_options": {"fo": path},
            "consolidated": False
        },
        chunks={}
    )
