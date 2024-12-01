import argparse
import os

import beaker
import beaker.services

from rollout import (
    EvaluatorRollout,
    ForcingDataset,
    InferenceRollout,
    InitialCondition,
    INITIAL_CONDITION_MOUNT_PATH,
    Model,
    ReferenceDataset,
)


IMAGE = "spencerc/fme-5190f642"
WANDB_RUN_GROUP = "ace2-som-ensemble"
MODELS = {
    "ACE2-SOM-multi-climate-RS0": "01J3R29H4B1XBMNC7WC2P624RF",
    "ACE2-SOM-multi-climate-RS1": "01J3YNZ3V1TZ2JB9YXYC1GM1CN",
    "ACE2-SOM-multi-climate-RS2": "01J484CQ8ZMRWDYG7ZD8DG77S1",
    "ACE2-SOM-multi-climate-RS3": "01J4BR6J5AW32ZDQ77VZ60P4KT",
    "ACE2-SOM-increasing-co2-RS0": "01JB5J1RJSBP3DXTPRD9N14AMR",
    "ACE2-SOM-increasing-co2-RS1": "01JB8NDKYPSSQ6TN7SMJR61DSZ",
    "ACE2-SOM-increasing-co2-RS2": "01JB8NHT9B5JY4DJ6A8M8741G0",
    "ACE2-SOM-increasing-co2-RS3": "01JBASP62H7K5NRZJ3GPY6V0DC"
}

# ACE2-SOM experiments
SPIN_UP_ROOT = "/climate-default/2024-08-15-vertically-resolved-1deg-c96-shield-som-ensemble-spin-up-fme-dataset/netcdfs"
FORCING_ROOT = "/climate-default/2024-07-09-vertically-resolved-1deg-c96-shield-som-ensemble-fme-dataset/netcdfs"
SPIN_UP_DATASETS = {
    "1xCO2": os.path.join(SPIN_UP_ROOT, "concatenated-1xCO2-ic_0005"),
    "2xCO2": os.path.join(SPIN_UP_ROOT, "concatenated-2xCO2-ic_0005"),
    "3xCO2": os.path.join(SPIN_UP_ROOT, "concatenated-3xCO2-ic_0002"),
    "4xCO2": os.path.join(SPIN_UP_ROOT, "concatenated-4xCO2-ic_0005"),
}
FORCING_DATASETS = {
    "1xCO2": os.path.join(FORCING_ROOT, "1xCO2-ic_0005"),
    "2xCO2": os.path.join(FORCING_ROOT, "2xCO2-ic_0005"),
    "3xCO2": os.path.join(FORCING_ROOT, "3xCO2-ic_0002"),
    "4xCO2": os.path.join(FORCING_ROOT, "4xCO2-ic_0005"),
}
INITIAL_CONDITIONS = {
    0: "2030-01-01T06:00:00",
    1: "2030-01-01T12:00:00",
    2: "2030-01-01T18:00:00",
    3: "2030-01-02T00:00:00",
    4: "2030-01-02T06:00:00",
}
ENSEMBLE_OVERLAY_TEMPLATE = {
    "n_forward_steps": 14600,
    "logging": {
        "log_to_screen": True,
        "log_to_wandb": True,
        "log_to_file": True,
        "project": "ace-som"
    },
    "forcing_loader": {
        "dataset": {"data_path": "/path/to/data"},
        "num_data_workers": 8
    },
    "initial_condition": {"path": "/initial_condition/restart.nc"},
    "data_writer": {
        "save_monthly_files": False,
        "save_prediction_files": False
    }
}

# SHiELD-SOM-C96 experiments
# Arbitrary one degree model with same outputs
ONE_DEGREE_MODEL = "01J4BR6J5AW32ZDQ77VZ60P4KT"
REFERENCE_DATASETS = {
    "1xCO2": os.path.join(FORCING_ROOT, "1xCO2-ic_0005"),
    "2xCO2": os.path.join(FORCING_ROOT, "2xCO2-ic_0005"),
    "3xCO2": os.path.join(FORCING_ROOT, "3xCO2-ic_0001"),
    "4xCO2": os.path.join(FORCING_ROOT, "4xCO2-ic_0005"),
}
REFERENCE_PREDICTION_ROOT = FORCING_ROOT

# SHiELD-SOM-C24 experiments
# Arbitrary four degree model with same outputs
FOUR_DEGREE_MODEL = "01J6BXZN8R1D3892JMWM7NRMDS"
FOUR_DEGREE_REFERENCE_ROOT = "/climate-default/2024-07-09-vertically-resolved-4deg-c96-shield-som-ensemble-fme-dataset/netcdfs"
FOUR_DEGREE_REFERENCE_DATASETS = {
    "1xCO2": os.path.join(FOUR_DEGREE_REFERENCE_ROOT, "1xCO2-ic_0005"),
    "2xCO2": os.path.join(FOUR_DEGREE_REFERENCE_ROOT, "2xCO2-ic_0005"),
    "3xCO2": os.path.join(FOUR_DEGREE_REFERENCE_ROOT, "3xCO2-ic_0001"),
    "4xCO2": os.path.join(FOUR_DEGREE_REFERENCE_ROOT, "4xCO2-ic_0005")
}
BASELINE_PREDICTION_ROOT = "/climate-default/2024-11-12-vertically-resolved-4deg-c24-shield-som-tuned-cdmbgwd-baseline-fme-dataset/netcdfs"
INITIAL_CONDITION_LABEL_MAPPING = {f"ic{n}": f"ic_{n:04d}" for n in range(1, 6)}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--spin-up", action="store_true")
    parser.add_argument("--experiments", action="store_true")
    args = parser.parse_args()

    client = beaker.Beaker.from_env()

    # Run the ACE2-SOM experiments
    for model_name, model_id in MODELS.items():
        model = Model(name=model_name, beaker_dataset=model_id)

        spin_up_experiments = {}
        for climate, spin_up_dataset in SPIN_UP_DATASETS.items():
            for offset, start_time in INITIAL_CONDITIONS.items():
                steps = 1460 - offset
                spin_up_initial_condition = InitialCondition(
                    name=f"ic{offset + 1}",
                    path=os.path.join(spin_up_dataset, "2030010100.nc"),
                    time=start_time,
                )
                spin_up_forcing = ForcingDataset(
                    name=climate, path=spin_up_dataset
                )
                spin_up_experiment = InferenceRollout(
                    client=client,
                    image=IMAGE,
                    model=model,
                    steps=steps,
                    forcing_dataset=spin_up_forcing,
                    initial_condition=spin_up_initial_condition,
                    tag="spin-up",
                    overlay_template=ENSEMBLE_OVERLAY_TEMPLATE
                )
                if args.spin_up:
                    spin_up_experiment.submit()
                spin_up_experiments[model_name, climate, offset] = spin_up_experiment

        if args.experiments:
            for (
                model_name,
                climate,
                offset,
            ), spin_up_experiment in spin_up_experiments.items():
                spun_up_ic = (
                    beaker.services.ExperimentClient(client)
                    .get(spin_up_experiment.name)
                    .jobs[-1]
                    .result.beaker
                )
                data_source = beaker.DataSource(beaker=spun_up_ic)
                initial_condition = InitialCondition(
                    name=f"ic{offset + 1}",
                    path=os.path.join(INITIAL_CONDITION_MOUNT_PATH, "restart.nc"),
                    time="2031-01-01T06:00:00",
                    data_source=data_source,
                )
                forcing = ForcingDataset(name=climate, path=FORCING_DATASETS[climate])
                model = Model(name=model_name, beaker_dataset=MODELS[model_name])
                rollout = InferenceRollout(
                    client=client,
                    image=IMAGE,
                    model=model,
                    steps=14604,
                    forcing_dataset=forcing,
                    initial_condition=initial_condition,
                    tag="ensemble-member",
                    overlay_template=ENSEMBLE_OVERLAY_TEMPLATE,
                    wandb_run_group=WANDB_RUN_GROUP,
                )
                rollout.submit()

    # Run the data-only experiments
    if args.experiments:
        model = Model(name="SHiELD-SOM-C96", beaker_dataset=ONE_DEGREE_MODEL)
        for (
            initial_condition_short,
            initial_condition_long,
        ) in INITIAL_CONDITION_LABEL_MAPPING.items():
            for climate, root in REFERENCE_DATASETS.items():
                steps = 14604
                reference_dataset = ReferenceDataset(
                    name=climate,
                    path=root,
                    start_times=["2031-01-01T06:00:00"],
                )
                prediction_dataset = ReferenceDataset(
                    name=f"{initial_condition_short}-{climate}",
                    path=os.path.join(
                        REFERENCE_PREDICTION_ROOT, f"{climate}-{initial_condition_long}"
                    ),
                    start_times=["2031-01-01T06:00:00"],
                )
                experiment = EvaluatorRollout(
                    client=client,
                    image=IMAGE,
                    model=model,
                    steps=steps,
                    reference_dataset=reference_dataset,
                    prediction_dataset=prediction_dataset,
                    wandb_run_group=WANDB_RUN_GROUP,
                )
                experiment.submit()

        model = Model(name="SHiELD-SOM-C24-tuned-cdmbgwd", beaker_dataset=FOUR_DEGREE_MODEL)
        for (
            initial_condition_short,
            initial_condition_long,
        ) in INITIAL_CONDITION_LABEL_MAPPING.items():
            for climate, root in FOUR_DEGREE_REFERENCE_DATASETS.items():
                steps = 14604
                reference_dataset = ReferenceDataset(
                    name=climate,
                    path=root,
                    start_times=["2031-01-01T06:00:00"],
                )
                prediction_dataset = ReferenceDataset(
                    name=f"{initial_condition_short}-{climate}",
                    path=os.path.join(
                        BASELINE_PREDICTION_ROOT, f"{climate}-{initial_condition_long}"
                    ),
                    start_times=["2031-01-01T06:00:00"],
                )
                experiment = EvaluatorRollout(
                    client=client,
                    image=IMAGE,
                    model=model,
                    steps=steps,
                    reference_dataset=reference_dataset,
                    prediction_dataset=prediction_dataset,
                    wandb_run_group=WANDB_RUN_GROUP,
                )
                experiment.submit()
