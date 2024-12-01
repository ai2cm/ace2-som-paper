# requires beaker-py, install with
# pip install -U beaker-py

import beaker
import copy
import uuid
from typing import Dict, Any
import tempfile
import yaml
import os


# Latest working commit of full-model; see GitHub issue #1286.
IMAGE_NAME = "spencerc/fme-29354918"  

LOCAL_BASE_CONFIG_FILENAME = "base-config-increasing-co2.yaml"
DATASET_CONFIG_FILENAME = "config.yaml"
DATASET_CONFIG_MOUNTPATH = "/configmount"

# Created the stats dataset somewhat manually, since our workflow in full-model
# doesn't work out of the box for excluding a period in the middle of a run.
STATS_DATASET_NAME = "spencerc/2024-10-25-vertically-resolved-1deg-c96-shield-som-stats-increasing-co2-combined-training-segments"
DATASET = "/climate-default/2024-07-16-vertically-resolved-1deg-c96-shield-som-increasing-co2-fme-dataset/netcdfs/increasing-CO2"


with open(LOCAL_BASE_CONFIG_FILENAME, "r") as file:
    BASE_CONFIG = yaml.safe_load(file)


NAMES = [
    "2024-10-26-ace-som-increasing-co2-RS0",
    "2024-10-26-ace-som-increasing-co2-RS1",
    "2024-10-26-ace-som-increasing-co2-RS2",
    "2024-10-26-ace-som-increasing-co2-RS3"
]


def write_config_dataset(config: Dict[str, Any]):
    with tempfile.TemporaryDirectory() as temp_dir:
        filepath = os.path.join(temp_dir, DATASET_CONFIG_FILENAME)
        with open(filepath, "w") as f:
            yaml.safe_dump(config, f)
        dataset_name = "ace-training-config-" + str(uuid.uuid4())[:8]
        dataset = client.dataset.create(dataset_name, filepath)
    return dataset


def get_experiment_spec(name: str, config: Dict[str, Any], image_name=IMAGE_NAME):
    """Given a dict representing the training configuration, return a beaker experiment spec."""
    config_dataset = write_config_dataset(config)
    env_vars = [
        beaker.EnvVar(name="WANDB_API_KEY", secret="wandb-api-key-ai2cm-sa"),
        beaker.EnvVar(name="WANDB_JOB_TYPE", value="training"),
        beaker.EnvVar(name="WANDB_NAME", value=name),
        beaker.EnvVar(name="NCCL_PROTO", value="^LL128"),
        beaker.EnvVar(name="WANDB_USERNAME", value="spencerc_ai2"),
        beaker.EnvVar(name="WANDB_RUN_GROUP", value="ace-som"),
    ]
    datasets = [
        beaker.DataMount(
            source=beaker.DataSource(beaker=config_dataset.id),
            mount_path=DATASET_CONFIG_MOUNTPATH,
        ),
        beaker.DataMount(
            source=beaker.DataSource(beaker=STATS_DATASET_NAME),
            mount_path="/statsdata",
        ),
        beaker.DataMount(
            mount_path="/climate-default",
            source=beaker.DataSource(weka="climate-default"),
        ),
    ]
    spec = beaker.ExperimentSpec(
        budget="ai2/climate",
        description="Train 1° ACE-SOM with a strict ACEv2-like configuration.",
        tasks=[
            beaker.TaskSpec(
                name=name,
                image=beaker.ImageSource(beaker=image_name),
                command=[
                    "torchrun",
                    "--nproc_per_node",
                    "8",
                    "-m",
                    "fme.ace.train",
                    f"{DATASET_CONFIG_MOUNTPATH}/{DATASET_CONFIG_FILENAME}",
                ],
                result=beaker.ResultSpec(path="/output"),
                resources=beaker.TaskResources(gpu_count=8, shared_memory="400GiB"),
                context=beaker.TaskContext(priority="high", preemptible=True),
                constraints=beaker.Constraints(cluster=["ai2/jupiter-cirrascale-2"]),
                env_vars=env_vars,
                datasets=datasets,
            )
        ],
    )
    return spec


if __name__ == "__main__":
    client = beaker.Beaker.from_env()
    for name in NAMES:
        config = copy.deepcopy(BASE_CONFIG)
        print(f"Creating experiment {name}.")
        spec = get_experiment_spec(name, config)
        experiment = client.experiment.create(name, spec, workspace="ai2/ace")
        print(f"Experiment created. See https://beaker.org/ex/{experiment.id}")
