import ujson

from kerchunk.netCDF3 import NetCDF3ToZarr


netcdf_file = snakemake.input.file
json_file = snakemake.output.json
chunks = NetCDF3ToZarr(netcdf_file)
with open(json_file, "wb") as file:
    output = ujson.dumps(chunks.translate()).encode()
    file.write(output)
