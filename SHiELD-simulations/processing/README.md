# Post-processing workflow for full model emulation data

This repository contains code for a `snakemake` workflow built to prepare
FV3GFS or SHiELD output for use in training and validating the Ai2 Climate
Emulator (ACE). It handles archiving the raw data, regridding raw output from
the native cubed sphere grid to the Gaussian grid, applying a spherical
harmonic transform round trip filter, rechunking, concatenating, converting to
zarr, and transferring to the cloud. `snakemake` allows the steps of this
workflow to be parallelized in SLURM batch jobs on the post-processing and data
transfer partitions of Gaea.

## Installation

The environment-related code has a soft assumption that one is working on
Gaea—namely it assumes that a system-installed `conda` is available and
installed in a particular hard-coded path. This could be relaxed to support
other HPC platforms. If on Gaea, however, installation is quite straightforward
through the `make create_environment` rule:

```
$ make create_environment
```

This creates an environment with the name assigned to the `ENVIRONMENT_NAME`
variable in the `Makefile`; this environment is automatically used when calling
other `make` rules.

## Usage

To apply the workflow to a run, or set of runs, modify the included
`config/config.yaml` file. This file contains the necessary metadata for the
workflow, including:

- Directories to archive data, write intermediate post-processing files, write
  final post-processed files, and upload data to the cloud.
- Root directories of simulations, assumed to be organized following typical
  SHiELD conventions, i.e. with `ascii`, `history`, and `restart`
  subdirectories, each with segment subdirectories within. These simulations
  need not have identical resolutions, segment lengths, total lengths, or
  outputs.
- A library of `fregrid`-compatible grid information for possible cubed sphere
  source grids used by the  simulations. 
- Target grid information, including paths to `fregrid`-compatible mosaic
  files, and grid-specific target chunk sizes for eventual zarr stores, since
  ideal chunk sizes depend on resolution.
- Configuration details for the when to apply the spherical harmonic transform
  round trip filter, including names of variables to exclude.

Once configured, it can be useful to run a "dry-run" of the workflow to see if
`snakemake` can find all the necessary input data for the workflow:

```
$ make dry_run_zarr
```

If this checks out, the full workflow can be run using:

```
$ make zarr transfer
```

This will execute things end-to-end, culminating with transferring the
generated zarr stores to the cloud. Note that it can be useful to submit this
via a `screen` session, as the workflow can take several hours to
complete—depending on the size of the zarr stores, the transfer process can
actually be the slowest part, since at the moment it is not particularly well
parallelized.

Note that we have leveraged the [job grouping
capability](https://snakemake.readthedocs.io/en/stable/executing/grouping.html)
of `snakemake` to combine groups of smaller jobs into longer running batch
jobs. This is used to address a common complaint of `snakemake`'s cluster
execution, which is that it generates an unhealthy number of batch jobs that
need to work through the queue. Still more could likely be done to optimize
this, however.
