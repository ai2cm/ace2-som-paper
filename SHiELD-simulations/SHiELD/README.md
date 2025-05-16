# `shield_run`

This repository contains code for configuring and running SHiELD simulations on
Gaea, through a similar pattern to what we use when running the Python wrapper
in the cloud.

Right now this repo mixes configuration and software, which we might consider
separating at some point—the code for running simulations we implemented is
general—but it is also meant to fully describe the process for configuring and
running the reference simulations we run for full model emulation with SHiELD,
for which this software is required.

## Installation

To install the necessary software, just a couple make commands are needed.
First we run `make build_shield` to compile SHiELD using GFDL's `SHiELD_build`
system. This automatically checks out any submodules as needed:

```
$ make build_shield
```

Then we create a Python environment using `make create_environment` within
which we install several Python packages used for configuring and running
simulations:

```
$ make create_environment
```

Note you can customize the name of the environment by editing the `Makefile`.

You can test that the software was installed properly by using the `make
test_shield_run` rule:

```
$ make test_shield_run
```

Note that this does not test running SHiELD—rather mainly helper code around
that—but we could consider adding some regression tests at some point.

## Running simulations

Running simulations is now straightforward.  The `python -m shield_run.submit`
tool is built for this purpose.  All that is required is the path to the SHiELD
executable, the path to the uncompiled configuration, the number of segments
the simulation will be run for, the wall clock time for each segment, and the
directory the simulation is to be completed within:

```
$ python -m shield_run.submit --help
usage: submit.py [-h] [--initial_conditions INITIAL_CONDITIONS]
                 [--fv3config_root FV3CONFIG_ROOT] [--resume]
                 executable config segments time destination

positional arguments:
  executable            Path to the SHiELD executable
  config                Path to the raw uncompiled config
  segments              Number of segments to run the simulation
  time                  Wall clock time per segment
  destination           Path to the output of the run

options:
  -h, --help            show this help message and exit
  --initial_conditions INITIAL_CONDITIONS
                        Path to the initial conditions for overriding what is
                        in the config (mainly for initial condition
                        ensembles).
  --fv3config_root FV3CONFIG_ROOT
                        Path to reference data for fv3config; used only for
                        appending grid files if requested.
  --resume              Whether to resume this simulation after an
                        interruption. Merely avoids the run creation step.
```

An example for submitting a run can be found in the `Makefile`.  The
configuration is defined in the `configs` directory.
