ENVS=$(shell pwd)/envs
CONDA_RUN=$(ENVS)/conda_run.sh

ENVIRONMENT_NAME=fme-post-processing
ENVIRONMENT_FILE=$(ENVS)/environment.yaml
UNCOMPILED_REQUIREMENTS=$(ENVS)/requirements.in
PIP_REQUIREMENTS=$(ENVS)/requirements.txt


create_environment:
	bash $(ENVS)/create_environment.sh \
		$(ENVIRONMENT_FILE) \
		$(ENVIRONMENT_NAME) \
		$(PIP_REQUIREMENTS)

lock_pip:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) pip-compile \
		--strip-extras \
		$(UNCOMPILED_REQUIREMENTS) > $(PIP_REQUIREMENTS)

dry_run_archive:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) snakemake \
		-n \
		--groups archive_diagnostics_file=archive_diagnostics \
			archive_scalar_diagnostics_file=archive_diagnostics \
			tar_restart_directory=archive_restarts \
			archive_restart_directory=archive_restarts \
			tar_ascii_directory=archive_ascii \
			archive_ascii_directory=archive_ascii \
		--group-components archive_diagnostics=288 archive_restarts=12 archive_ascii=12 \
		--executor cluster-generic \
		--cluster-generic-submit-cmd "sbatch --mem={resources.mem} --time=16:00:00 --partition={resources.partition} --constraint={resources.constraint} --output=slurm-logs/%x.%j.out" \
		--cluster-generic-cancel-cmd "scancel" \
		--jobs 75 \
		archive

archive:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) snakemake \
		--groups archive_diagnostics_file=archive_diagnostics \
			archive_scalar_diagnostics_file=archive_diagnostics \
			tar_restart_directory=archive_restarts \
			archive_restart_directory=archive_restarts \
			tar_ascii_directory=archive_ascii \
			archive_ascii_directory=archive_ascii \
		--group-components archive_diagnostics=144 archive_restarts=12 archive_ascii=12 \
		--executor cluster-generic \
		--cluster-generic-submit-cmd "sbatch --mem={resources.mem} --time=16:00:00 --partition={resources.partition} --constraint={resources.constraint} --output=slurm-logs/%x.%j.out" \
		--cluster-generic-cancel-cmd "scancel" \
		--jobs 75 \
		archive

dry_run_zarr:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) snakemake \
		-n \
		--groups regrid_tape=regrid kerchunk_regridded_segment_tape=kerchunk kerchunk_scalar_segment_tape=kerchunk_scalar_tape \
		--group-components regrid=10 kerchunk=100 kerchunk_scalar_tape=100 \
		--executor cluster-generic \
		--cluster-generic-submit-cmd "sbatch --mem={resources.mem} --time=10:00:00 --partition={resources.partition} --constraint={resources.constraint} --output=slurm-logs/%x.%j.out" \
		--cluster-generic-cancel-cmd "scancel" \
		--jobs 75 \
		zarr

zarr:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) snakemake \
		--groups regrid_tape=regrid kerchunk_regridded_segment_tape=kerchunk kerchunk_scalar_segment_tape=kerchunk_scalar_tape \
		--group-components regrid=10 kerchunk=100 kerchunk_scalar_tape=100 \
		--executor cluster-generic \
		--cluster-generic-submit-cmd "sbatch --mem={resources.mem} --time=16:00:00 --partition={resources.partition} --constraint={resources.constraint} --output=slurm-logs/%x.%j.out" \
		--cluster-generic-cancel-cmd "scancel" \
		--jobs 75 \
		zarr

transfer:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) snakemake \
		--executor cluster-generic \
		--cluster-generic-submit-cmd "sbatch --mem={resources.mem} --time=16:00:00 --partition={resources.partition} --constraint={resources.constraint} --output=slurm-logs/%x.%j.out" \
		--cluster-generic-cancel-cmd "scancel" \
		--jobs 75 \
		transfer


unlock:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) snakemake --unlock
