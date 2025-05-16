# Variables one might modify
SCRATCH_ROOT=/gpfs/f5/gfdl_w/scratch/Spencer.Clark/scratch/reference-SHiELD-runs
INITIAL_CONDITIONS_ROOT=$(SCRATCH_ROOT)/initial-conditions
GRID_DATA_ROOT=/gpfs/f5/gfdl_w/scratch/Spencer.Clark/reference/SHiELD-SOM-input-data/input-data/data/grid_data/v1.2
ENVIRONMENT_NAME=SHiELD-simulations
CONFIG_DIR=$(shell pwd)/config

# Static variables
SOFTWARE_DIR=$(shell pwd)/software
SCRIPTS_DIR=$(shell pwd)/scripts
ENVIRONMENT_FILE=$(SOFTWARE_DIR)/environment.yaml
PIP_REQUIREMENTS=$(SOFTWARE_DIR)/requirements.txt
CONDA_RUN=$(SOFTWARE_DIR)/conda_run.sh
EXECUTABLE=$(SOFTWARE_DIR)/submodules/SHiELD_build/Build/bin/SHiELD_nh.repro.64bit.intel.x

# A Python package for running SHiELD simulations via fv3config and SLURM
SHIELD_RUN_DIR=$(SOFTWARE_DIR)/shield_run
SHIELD_RUN_TESTS=$(SHIELD_RUN_DIR)/shield_run/tests

# Makefile-specific variables
CLIM_SST_INITIAL_CONDITIONS=0001 0002
SOM_INITIAL_CONDITIONS=0001 0002 0003 0004 0005
CLIMATES=1xCO2 2xCO2 3xCO2 4xCO2

# Key-value pairs to infer initial condition date from initial condition name
# https://stackoverflow.com/questions/62005888/key-value-pair-in-makefile
IC_DATE_CLIM_SST_0001 := 2034093000
IC_DATE_CLIM_SST_0002 := 2034100100

IC_DATE_C96_ic_0001 := 2029122800
IC_DATE_C96_ic_0002 := 2029122900
IC_DATE_C96_ic_0003 := 2029123000
IC_DATE_C96_ic_0004 := 2029123100
IC_DATE_C96_ic_0005 := 2030010100

IC_DATE_C96_3xCO2_ic_0001 := 2029123100
IC_DATE_C96_3xCO2_ic_0002 := 2030010100
IC_DATE_C96_3xCO2_ic_0003 := 2029123000
IC_DATE_C96_3xCO2_ic_0004 := 2029122900
IC_DATE_C96_3xCO2_ic_0005 := 2029122800

IC_DATE_C24_0001 := 2030010100
IC_DATE_C24_0002 := 2029123100
IC_DATE_C24_0003 := 2029123000
IC_DATE_C24_0004 := 2029122900
IC_DATE_C24_0005 := 2029122800

update_submodules:
	git submodule update --init --recursive

build_shield:
	bash $(SOFTWARE_DIR)/build_SHiELD.sh

create_environment:
	bash $(SOFTWARE_DIR)/create_environment.sh \
		$(ENVIRONMENT_FILE) \
		$(ENVIRONMENT_NAME) \
		$(PIP_REQUIREMENTS) \
		$(SHIELD_RUN_DIR)

lock_pip:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) pip-compile \
		$(ENVIRONMENT_DIR)/requirements.in > $(PIP_REQUIREMENTS)

reformat:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) black $(SHIELD_RUN_DIR)

test_shield_run:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) pytest -vv $(SHIELD_RUN_TESTS)

submit_repeating_sst_spin_up_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/repeating-sst-spin-up-C24-tuned-cdmbgwd.yaml \
			44 \
			01:20:00 \
			$(SCRATCH_ROOT)/repeating-sst-spin-up-C24-tuned-cdmbgwd \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_repeating_sst_bridge_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/repeating-sst-bridge-C24-tuned-cdmbgwd.yaml \
			1 \
			01:20:00 \
			$(SCRATCH_ROOT)/repeating-sst-bridge-C24-tuned-cdmbgwd \
			--initial-conditions=$(SCRATCH_ROOT)/repeating-sst-spin-up-C24-tuned-cdmbgwd/restart/2034090100 \
			--grid-data-root=$(GRID_DATA_ROOT)

prepare_repeating_sst_C24_tuned_cdmbgwd_initial_conditions:
	$(SCRIPTS_DIR)/prepare_restart_files.sh \
		$(SCRATCH_ROOT)/repeating-sst-bridge-C24-tuned-cdmbgwd/restart/2034100100 \
		$(INITIAL_CONDITIONS_ROOT)/repeating-sst-bridge-C24-tuned-cdmbgwd \
		2034093000 \
		2034100100

submit_repeating_sst_C24_tuned_cdmbgwd: $(addprefix submit_repeating_sst_C24_tuned_cdmbgwd_ic_, $(CLIM_SST_INITIAL_CONDITIONS))
submit_repeating_sst_C24_tuned_cdmbgwd_ic_%:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/repeating-sst-C24-tuned-cdmbgwd.yaml \
			61 \
			01:20:00 \
			$(SCRATCH_ROOT)/repeating-sst-C24-tuned-cdmbgwd-ic-$* \
			--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/repeating-sst-bridge-C24-tuned-cdmbgwd/$(IC_DATE_CLIM_SST_$*) \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_repeating_sst_spin_up_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/repeating-sst-spin-up-C96.yaml \
			177 \
			00:42:00 \
			$(SCRATCH_ROOT)/repeating-sst-spin-up-C96 \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_repeating_sst_bridge_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/repeating-sst-bridge-C96.yaml \
			1 \
			00:42:00 \
			$(SCRATCH_ROOT)/repeating-sst-bridge-C96 \
			--initial-conditions=$(SCRATCH_ROOT)/repeating-sst-spin-up-C96/restart/2034090100 \
			--grid-data-root=$(GRID_DATA_ROOT)

prepare_repeating_sst_C96_initial_conditions:
	$(SCRIPTS_DIR)/prepare_restart_files.sh \
		$(SCRATCH_ROOT)/repeating-sst-bridge-C96/restart/2034100100 \
		$(INITIAL_CONDITIONS_ROOT)/repeating-sst-bridge-C96 \
		2034093000 \
		2034100100

submit_repeating_sst_C96: $(addprefix submit_repeating_sst_C96_, $(CLIM_SST_INITIAL_CONDITIONS))
submit_repeating_sst_C96_%:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/repeating-sst-C96.yaml \
			183 \
			00:42:00 \
			$(SCRATCH_ROOT)/repeating-sst-C96-$* \
			--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/repeating-sst-bridge-C96/$(IC_DATE_CLIM_SST_$*)
			--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_spin_up_1xCO2_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-spin-up-1xCO2-C24-tuned-cdmbgwd.yaml \
			30 \
			01:20:00 \
			$(SCRATCH_ROOT)/slab-ocean-spin-up-1xCO2-C24-tuned-cdmbgwd \
			--initial-conditions=$(SCRATCH_ROOT)/repeating-sst-C24-tuned-cdmbgwd-ic_0002/restart/2050010100 \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_bridge_1xCO2_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-bridge-1xCO2-C24-tuned-cdmbgwd.yaml \
			1 \
			01:20:00 \
			$(SCRATCH_ROOT)/slab-ocean-bridge-1xCO2-C24-tuned-cdmbgwd \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-1xCO2-C24-tuned-cdmbgwd/restart/2029090100 \
			--grid-data-root=$(GRID_DATA_ROOT)

prepare_slab_ocean_1xCO2_C24_tuned_cdmbgwd_initial_conditions:
	$(SCRIPTS_DIR)/prepare_restart_files.sh \
		$(SCRATCH_ROOT)/slab-ocean-bridge-1xCO2-C24-tuned-cdmbgwd/restart/2030010100 \
		$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-1xCO2-C24-tuned-cdmbgwd \
		2029122600 \
		2029122700 \
		2029122800 \
		2029122900 \
		2029123000 \
		2029123100 \
		2030010100

submit_slab_ocean_spin_up_2xCO2_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-spin-up-2xCO2-C24-tuned-cdmbgwd.yaml \
			30 \
			01:20:00 \
			$(SCRATCH_ROOT)/slab-ocean-spin-up-2xCO2-C24-tuned-cdmbgwd \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-1xCO2-C24-tuned-cdmbgwd/restart/2030010100 \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_bridge_2xCO2_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-bridge-2xCO2-C24-tuned-cdmbgwd.yaml \
			1 \
			01:20:00 \
			$(SCRATCH_ROOT)/slab-ocean-bridge-2xCO2-C24-tuned-cdmbgwd \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-2xCO2-C24-tuned-cdmbgwd/restart/2029090100 \
			--grid-data-root=$(GRID_DATA_ROOT)

prepare_slab_ocean_2xCO2_C24_tuned_cdmbgwd_initial_conditions:
	$(SCRIPTS_DIR)/prepare_restart_files.sh \
		$(SCRATCH_ROOT)/slab-ocean-bridge-2xCO2-C24-tuned-cdmbgwd/restart/2030010100 \
		$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-2xCO2-C24-tuned-cdmbgwd \
		2029122800 \
		2029122900 \
		2029123000 \
		2029123100 \
		2030010100

submit_slab_ocean_spin_up_3xCO2_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-spin-up-3xCO2-C24-tuned-cdmbgwd.yaml \
			30 \
			01:20:00 \
			$(SCRATCH_ROOT)/slab-ocean-spin-up-3xCO2-C24-tuned-cdmbgwd \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-1xCO2-C24-tuned-cdmbgwd/restart/2030010100 \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_bridge_3xCO2_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-bridge-3xCO2-C24-tuned-cdmbgwd.yaml \
			1 \
			01:20:00 \
			$(SCRATCH_ROOT)/slab-ocean-bridge-3xCO2-C24-tuned-cdmbgwd \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-3xCO2-C24-tuned-cdmbgwd/restart/2029090100 \
			--grid-data-root=$(GRID_DATA_ROOT)

prepare_slab_ocean_3xCO2_C24_tuned_cdmbgwd_initial_conditions:
	$(SCRIPTS_DIR)/prepare_restart_files.sh \
		$(SCRATCH_ROOT)/slab-ocean-bridge-3xCO2-C24-tuned-cdmbgwd/restart/2030010100 \
		$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-3xCO2-C24-tuned-cdmbgwd \
		2029122800 \
		2029122900 \
		2029123000 \
		2029123100 \
		2030010100

submit_slab_ocean_spin_up_4xCO2_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-spin-up-4xCO2-C24-tuned-cdmbgwd.yaml \
			30 \
			01:20:00 \
			$(SCRATCH_ROOT)/slab-ocean-spin-up-4xCO2-C24-tuned-cdmbgwd \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-1xCO2-C24-tuned-cdmbgwd/restart/2030010100 \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_bridge_4xCO2_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-bridge-4xCO2-C24-tuned-cdmbgwd.yaml \
			1 \
			01:20:00 \
			$(SCRATCH_ROOT)/slab-ocean-bridge-4xCO2-C24-tuned-cdmbgwd \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-4xCO2-C24-tuned-cdmbgwd/restart/2029090100 \
			--grid-data-root=$(GRID_DATA_ROOT)

prepare_slab_ocean_4xCO2_C24_tuned_cdmbgwd_initial_conditions:
	$(SCRIPTS_DIR)/prepare_restart_files.sh \
		$(SCRATCH_ROOT)/slab-ocean-bridge-4xCO2-C24-tuned-cdmbgwd/restart/2030010100 \
		$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-4xCO2-C24-tuned-cdmbgwd \
		2029122800 \
		2029122900 \
		2029123000 \
		2029123100 \
		2030010100

submit_SOM_production_1xCO2_C24_tuned_cdmbgwd: $(addprefix submit_SOM_production_1xCO2_C24_tuned_cdmbgwd_ic_, $(SOM_INITIAL_CONDITIONS))
submit_SOM_production_1xCO2_C24_tuned_cdmbgwd_ic_%:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
		$(EXECUTABLE) \
		$(CONFIG_DIR)/slab-ocean-1xCO2-C24-tuned-cdmbgwd.yaml \
		33 \
		01:20:00 \
		$(SCRATCH_ROOT)/slab-ocean-1xCO2-C24-tuned-cdmbgwd-ic_$* \
		--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-1xCO2-C24-tuned-cdmbgwd/$(IC_DATE_C24_$*) \
		--grid-data-root=$(GRID_DATA_ROOT)

submit_SOM_production_1xCO2_C24_tuned_cdmbgwd_radiation_multi_call:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
		$(EXECUTABLE) \
		$(CONFIG_DIR)/slab-ocean-1xCO2-C24-tuned-cdmbgwd-radiation-multi-call.yaml \
		6 \
		01:45:00 \
		$(SCRATCH_ROOT)/slab-ocean-1xCO2-C24-tuned-cdmbgwd-radiation-multi-call \
		--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-1xCO2-C24-tuned-cdmbgwd/2029122600 \
		--grid-data-root=$(GRID_DATA_ROOT)

submit_SOM_production_2xCO2_C24_tuned_cdmbgwd: $(addprefix submit_SOM_production_2xCO2_C24_tuned_cdmbgwd_ic_, $(SOM_INITIAL_CONDITIONS))
submit_SOM_production_2xCO2_C24_tuned_cdmbgwd_ic_%:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
		$(EXECUTABLE) \
		$(CONFIG_DIR)/slab-ocean-2xCO2-C24-tuned-cdmbgwd.yaml \
		33 \
		01:20:00 \
		$(SCRATCH_ROOT)/slab-ocean-2xCO2-C24-tuned-cdmbgwd-ic_$* \
		--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-2xCO2-C24-tuned-cdmbgwd/$(IC_DATE_C24_$*) \
		--grid-data-root=$(GRID_DATA_ROOT)

submit_SOM_production_3xCO2_C24_tuned_cdmbgwd: $(addprefix submit_SOM_production_3xCO2_C24_tuned_cdmbgwd_ic_, $(SOM_INITIAL_CONDITIONS))
submit_SOM_production_3xCO2_C24_tuned_cdmbgwd_ic_%:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
		$(EXECUTABLE) \
		$(CONFIG_DIR)/slab-ocean-3xCO2-C24-tuned-cdmbgwd.yaml \
		33 \
		01:20:00 \
		$(SCRATCH_ROOT)/slab-ocean-3xCO2-C24-tuned-cdmbgwd-ic_$* \
		--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-3xCO2-C24-tuned-cdmbgwd/$(IC_DATE_C24_$*) \
		--grid-data-root=$(GRID_DATA_ROOT)

submit_SOM_production_4xCO2_C24_tuned_cdmbgwd: $(addprefix submit_SOM_production_4xCO2_C24_tuned_cdmbgwd_ic_, $(SOM_INITIAL_CONDITIONS))
submit_SOM_production_4xCO2_C24_tuned_cdmbgwd_ic_%:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
		$(EXECUTABLE) \
		$(CONFIG_DIR)/slab-ocean-4xCO2-C24-tuned-cdmbgwd.yaml \
		33 \
		01:20:00 \
		$(SCRATCH_ROOT)/slab-ocean-4xCO2-C24-tuned-cdmbgwd-ic_$* \
		--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-4xCO2-C24-tuned-cdmbgwd/$(IC_DATE_C24_$*) \
		--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_spin_up_1xCO2_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-spin-up-1xCO2-C96.yaml \
			120 \
			00:42:00 \
			$(SCRATCH_ROOT)/slab-ocean-spin-up-1xCO2-C96 \
			--initial-conditions=$(SCRATCH_ROOT)/repeating-sst-C96-ic_0002/restart/2050010100 \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_bridge_1xCO2_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-bridge-1xCO2-C96.yaml \
			1 \
			00:42:00 \
			$(SCRATCH_ROOT)/slab-ocean-bridge-1xCO2-C96 \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-1xCO2-C96/restart/2029120100 \
			--grid-data-root=$(GRID_DATA_ROOT)

prepare_slab_ocean_1xCO2_C96_initial_conditions:
	$(SCRIPTS_DIR)/prepare_restart_files.sh \
		$(SCRATCH_ROOT)/slab-ocean-bridge-1xCO2-C96/restart/2030010100 \
		$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-1xCO2-C96 \
		2029122600 \
		2029122700 \
		2029122800 \
		2029122900 \
		2029123000 \
		2029123100 \
		2030010100

submit_slab_ocean_spin_up_2xCO2_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-spin-up-2xCO2-C96.yaml \
			120 \
			00:42:00 \
			$(SCRATCH_ROOT)/slab-ocean-spin-up-2xCO2-C24 \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-1xCO2-C96/restart/2030010100 \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_bridge_2xCO2_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-bridge-2xCO2-C96.yaml \
			1 \
			00:42:00 \
			$(SCRATCH_ROOT)/slab-ocean-bridge-2xCO2-C96 \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-2xCO2-C96/restart/2029120100 \
			--grid-data-root=$(GRID_DATA_ROOT)

prepare_slab_ocean_2xCO2_C96_initial_conditions:
	$(SCRIPTS_DIR)/prepare_restart_files.sh \
		$(SCRATCH_ROOT)/slab-ocean-bridge-2xCO2-C96/restart/2030010100 \
		$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-2xCO2-C96 \
		2029122800 \
		2029122900 \
		2029123000 \
		2029123100 \
		2030010100

submit_slab_ocean_spin_up_3xCO2_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-spin-up-3xCO2-C96.yaml \
			120 \
			00:42:00 \
			$(SCRATCH_ROOT)/slab-ocean-spin-up-3xCO2-C96 \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-1xCO2-C96/restart/2030010100 \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_bridge_3xCO2_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-bridge-3xCO2-C96.yaml \
			1 \
			00:42:00 \
			$(SCRATCH_ROOT)/slab-ocean-bridge-3xCO2-C96 \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-3xCO2-C96/restart/2029120100 \
			--grid-data-root=$(GRID_DATA_ROOT)

prepare_slab_ocean_3xCO2_C96_initial_conditions:
	$(SCRIPTS_DIR)/prepare_restart_files.sh \
		$(SCRATCH_ROOT)/slab-ocean-bridge-3xCO2-C96/restart/2030010100 \
		$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-3xCO2-C96 \
		2029122800 \
		2029122900 \
		2029123000 \
		2029123100 \
		2030010100

submit_slab_ocean_spin_up_4xCO2_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-spin-up-4xCO2-C96.yaml \
			120 \
			00:42:00 \
			$(SCRATCH_ROOT)/slab-ocean-spin-up-4xCO2-C96 \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-1xCO2-C96/restart/2030010100 \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_bridge_4xCO2_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-bridge-4xCO2-C96.yaml \
			1 \
			00:42:00 \
			$(SCRATCH_ROOT)/slab-ocean-bridge-4xCO2-C96 \
			--initial-conditions=$(SCRATCH_ROOT)/slab-ocean-spin-up-4xCO2-C96/restart/2029120100 \
			--grid-data-root=$(GRID_DATA_ROOT)

prepare_slab_ocean_4xCO2_C96_initial_conditions:
	$(SCRIPTS_DIR)/prepare_restart_files.sh \
		$(SCRATCH_ROOT)/slab-ocean-bridge-4xCO2-C96/restart/2030010100 \
		$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-4xCO2-C96 \
		2029122800 \
		2029122900 \
		2029123000 \
		2029123100 \
		2030010100

submit_SOM_production_1xCO2_C96: $(addprefix submit_SOM_production_1xCO2_C96_, $(SOM_INITIAL_CONDITIONS))
submit_SOM_production_1xCO2_C96_%:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
		$(EXECUTABLE) \
		$(CONFIG_DIR)/slab-ocean-1xCO2-C96.yaml \
		132 \
		00:42:00 \
		$(SCRATCH_ROOT)/slab-ocean-1xCO2-C96-$* \
		--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-1xCO2-C96/$(IC_DATE_C96_$*) \
		--grid-data-root=$(GRID_DATA_ROOT)

submit_SOM_production_1xCO2_C96_radiation_multi_call:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
		$(EXECUTABLE) \
		$(CONFIG_DIR)/slab-ocean-1xCO2-C96-radiation-multi-call.yaml \
		24 \
		01:30:00 \
		$(SCRATCH_ROOT)/slab-ocean-1xCO2-C96-radiation-multi-call \
		--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-1xCO2-C96/2029122600 \
		--grid-data-root=$(GRID_DATA_ROOT)

submit_SOM_production_2xCO2_C96: $(addprefix submit_SOM_production_2xCO2_C96_, $(SOM_INITIAL_CONDITIONS))
submit_SOM_production_2xCO2_C96_%:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
		$(EXECUTABLE) \
		$(CONFIG_DIR)/slab-ocean-2xCO2-C96.yaml \
		132 \
		00:42:00 \
		$(SCRATCH_ROOT)/slab-ocean-2xCO2-C96-$* \
		--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-2xCO2-C96/$(IC_DATE_C96_$*) \
		--grid-data-root=$(GRID_DATA_ROOT)

submit_SOM_production_3xCO2_C96: $(addprefix submit_SOM_production_3xCO2_C96_, $(SOM_INITIAL_CONDITIONS))
submit_SOM_production_3xCO2_C96_%:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
		$(EXECUTABLE) \
		$(CONFIG_DIR)/slab-ocean-3xCO2-C96.yaml \
		132 \
		00:42:00 \
		$(SCRATCH_ROOT)/slab-ocean-3xCO2-C96-$* \
		--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-3xCO2-C96/$(IC_DATE_C96_3xCO2_$*) \
		--grid-data-root=$(GRID_DATA_ROOT)

submit_SOM_production_4xCO2_C96: $(addprefix submit_SOM_production_4xCO2_C96_, $(SOM_INITIAL_CONDITIONS))
submit_SOM_production_4xCO2_C96_%:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
		$(EXECUTABLE) \
		$(CONFIG_DIR)/slab-ocean-4xCO2-C96.yaml \
		132 \
		00:42:00 \
		$(SCRATCH_ROOT)/slab-ocean-4xCO2-C96-$* \
		--initial-conditions=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-4xCO2-C96/$(IC_DATE_C96_$*) \
		--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_increasing_CO2_C24_tuned_cdmbgwd:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-increasing-CO2-C24-tuned-cdmbgwd.yaml \
			213 \
			01:20:00 \
			$(SCRATCH_ROOT)/slab-ocean-increasing-CO2-C24-tuned-cdmbgwd \
			--initial-condition=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-1xCO2-C24-tuned-cdmbgwd/2029122700 \
			--grid-data-root=$(GRID_DATA_ROOT)

submit_slab_ocean_increasing_CO2_C96:
	$(CONDA_RUN) $(ENVIRONMENT_NAME) python -m shield_run.submit \
			$(EXECUTABLE) \
			$(CONFIG_DIR)/slab-ocean-increasing-CO2-C96.yaml \
			852 \
			00:42:00 \
			$(SCRATCH_ROOT)/slab-ocean-increasing-CO2-C96 \
			--initial-condition=$(INITIAL_CONDITIONS_ROOT)/slab-ocean-bridge-1xCO2-C96/2029122700 \
			--grid-data-root=$(GRID_DATA_ROOT)
