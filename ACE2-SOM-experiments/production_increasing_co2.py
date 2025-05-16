import beaker

from rollout import (
    EvaluatorRollout,
    Model,
    ReferenceDataset,
)

IMAGE = "spencerc/fme-5190f642"
WANDB_RUN_GROUP = "ace2-som-increasing-co2"
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
ONE_DEGREE_MODEL = "01J4BR6J5AW32ZDQ77VZ60P4KT"
FOUR_DEGREE_MODEL = "01J6BXZN8R1D3892JMWM7NRMDS"

REFERENCE_DATASET = "/climate-default/2024-07-16-vertically-resolved-1deg-c96-shield-som-increasing-co2-fme-dataset/netcdfs/increasing-CO2"
FOUR_DEGREE_REFERENCE_DATASET = "/climate-default/2024-07-16-vertically-resolved-4deg-c96-shield-som-increasing-co2-fme-dataset/netcdfs/increasing-CO2"
BASELINE_DATASET = "/climate-default/2024-11-12-vertically-resolved-4deg-c24-shield-som-tuned-cdmbgwd-baseline-fme-dataset/netcdfs/increasing-CO2"

EVALUATOR_OVERLAY_TEMPLATE = {
    "n_forward_steps": 102267,
    "loader": {
        "start_indices": {"times": ["2031-01-01T06:00:00"]},
        "dataset": {"data_path": REFERENCE_DATASET},
        "num_data_workers": 4,
    },
    "data_writer": {"save_prediction_files": False, "save_monthly_files": False},
    "aggregator": {"log_zonal_mean_images": False, "log_histograms": True},
}


if __name__ == "__main__":
    client = beaker.Beaker.from_env()

    steps = 102267
    for model_name, model_id in MODELS.items():
        model = Model(name=model_name, beaker_dataset=model_id)
        reference_dataset = ReferenceDataset(
            name="increasing-co2",
            path=REFERENCE_DATASET,
            start_times=["2031-01-01T06:00:00"],
        )
        experiment = EvaluatorRollout(
            client=client,
            image=IMAGE,
            model=model,
            steps=steps,
            reference_dataset=reference_dataset,
            overlay_template=EVALUATOR_OVERLAY_TEMPLATE,
            wandb_run_group=WANDB_RUN_GROUP,
            priority="high"
        )
        experiment.submit()

    model = Model(name="SHiELD-SOM-C96", beaker_dataset=ONE_DEGREE_MODEL)
    reference_dataset = ReferenceDataset(
        name="increasing-co2",
        path=REFERENCE_DATASET,
        start_times=["2031-01-01T06:00:00"],
    )
    prediction_dataset = ReferenceDataset(
        name="increasing-co2",
        path=REFERENCE_DATASET,
        start_times=["2031-01-01T06:00:00"],
    )
    experiment = EvaluatorRollout(
        client=client,
        image=IMAGE,
        model=model,
        steps=steps,
        reference_dataset=reference_dataset,
        prediction_dataset=prediction_dataset,
        overlay_template=EVALUATOR_OVERLAY_TEMPLATE,
        wandb_run_group=WANDB_RUN_GROUP,
        priority="high"
    )
    experiment.submit()

    model = Model(name="SHiELD-SOM-C24-tuned-cdmbgwd", beaker_dataset=FOUR_DEGREE_MODEL)
    reference_dataset = ReferenceDataset(
        name="increasing-co2",
        path=FOUR_DEGREE_REFERENCE_DATASET,
        start_times=["2031-01-01T06:00:00"],
    )
    prediction_dataset = ReferenceDataset(
        name="increasing-co2",
        path=BASELINE_DATASET,
        start_times=["2031-01-01T06:00:00"],
    )
    experiment = EvaluatorRollout(
        client=client,
        image=IMAGE,
        model=model,
        steps=steps,
        reference_dataset=reference_dataset,
        prediction_dataset=prediction_dataset,
        overlay_template=EVALUATOR_OVERLAY_TEMPLATE,
        wandb_run_group=WANDB_RUN_GROUP,
        priority="high"
    )
    experiment.submit()
