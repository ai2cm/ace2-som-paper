# requires beaker-py, install with
# pip install -U beaker-py

import beaker
import copy
import uuid
from typing import Dict, Any
import tempfile
import yaml
import os


IMAGE_NAME = "spencerc/fme-af79580d"
LOCAL_BASE_CONFIG_FILENAME = "base-config-multi-climate.yaml"
DATASET_CONFIG_FILENAME = "config.yaml"
DATASET_CONFIG_MOUNTPATH = "/configmount"

# Note due to full-model#944, the description of the Beaker dataset
# associated with this stats dataset makes it look as though data
# from the 3xCO2 was included, but it was not.  This dataset was 
# constructed only from ensemble members 0001 through 0004 of the 
# 1xCO2, 2xCO2, and 4xCO2 climates.
STATS_DATASET_NAME = "andrep/2024-07-09-vertically-resolved-1deg-fme-c96-shield-som-ensemble-dataset-stats"


with open(LOCAL_BASE_CONFIG_FILENAME, "r") as file:
    BASE_CONFIG = yaml.safe_load(file)


NAMES = [
    "2024-07-23-ace-som-RS0-1deg-multi-climate-inline-inference",
    "2024-07-23-ace-som-RS1-1deg-multi-climate-inline-inference",
    "2024-07-23-ace-som-RS2-1deg-multi-climate-inline-inference",
    "2024-07-23-ace-som-RS3-1deg-multi-climate-inline-inference",
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
