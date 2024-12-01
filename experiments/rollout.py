import copy
import dataclasses
import datetime
import os
import tempfile
import uuid

import beaker
import beaker.services
import yaml

from typing import Dict, List, Optional

CLUSTERS = ["ai2/jupiter-cirrascale-2"]

WANDB_JOB_TYPE = "inference"
WANDB_API_KEY = "wandb-api-key-ai2cm-sa"
WANDB_USERNAME = "spencerc_ai2"
WANDB_PROJECT = "ace-som"

EVALUATOR = "fme.ace.evaluator"
INFERENCE = "fme.ace.inference"

CHECKPOINT_NAME = "best_inference_ckpt.tar"
BASE_CONFIG_NAME = "base-config-minimal.yaml"
DATASET_CONFIG_FILENAME = "config.yaml"
DATASET_CONFIG_MOUNTPATH = "/configmount"
INITIAL_CONDITION_MOUNT_PATH = "/initial_condition"

DEFAULT_OVERLAY_TEMPLATE = {
    "n_forward_steps": 14604,
    "logging": {"log_to_screen": True, "log_to_wandb": True, "log_to_file": True, "project": "ace-som"},
    "loader": {
        "start_indices": {"times": ["2020-01-01T06:00:00"]},
        "dataset": {"data_path": "/path/to/data"},
        "num_data_workers": 8,
    },
    "data_writer": {"save_monthly_files": False, "save_prediction_files": False},
}


@dataclasses.dataclass
class ReferenceDataset:
    name: str
    path: str
    start_times: List[str]


@dataclasses.dataclass
class Model:
    name: str
    beaker_dataset: str


class BaseRollout:
    @property
    def name(self):
        raise NotImplementedError

    @property
    def overlay(self):
        raise NotImplementedError

    @property
    def datasets(self):
        raise NotImplementedError

    @property
    def entrypoint(self):
        raise NotImplementedError

    @property
    def config(self):
        with open(BASE_CONFIG_NAME, "r") as f:
            base_config = yaml.safe_load(f)
        return {**base_config, **self.overlay}

    def write_config_dataset(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            filepath = os.path.join(temp_dir, DATASET_CONFIG_FILENAME)
            with open(filepath, "w") as f:
                yaml.safe_dump(self.config, f)
            dataset_name = "ace-inference-config-" + str(uuid.uuid4())[:8]
            dataset = self.client.dataset.create(dataset_name, filepath)
        return dataset

    @property
    def standard_datasets(self):
        config_dataset = self.write_config_dataset()
        datasets = [
            beaker.DataMount(
                source=beaker.DataSource(beaker=config_dataset.id),
                mount_path=DATASET_CONFIG_MOUNTPATH,
            ),
            beaker.DataMount(
                mount_path="/ckpt.tar",
                source=beaker.DataSource(beaker=self.model.beaker_dataset),
                sub_path=f"training_checkpoints/{CHECKPOINT_NAME}",
            ),
            beaker.DataMount(
                mount_path="/climate-default",
                source=beaker.DataSource(weka="climate-default"),
            ),
        ]
        return datasets

    @property
    def experiment_spec(self):
        env_vars = [
            beaker.EnvVar(name="WANDB_API_KEY", secret=WANDB_API_KEY),
            beaker.EnvVar(name="WANDB_JOB_TYPE", value=WANDB_JOB_TYPE),
            beaker.EnvVar(name="WANDB_USERNAME", value=WANDB_USERNAME),
            beaker.EnvVar(name="WANDB_PROJECT", value=WANDB_PROJECT),
            beaker.EnvVar(name="WANDB_NAME", value=self.name),
            beaker.EnvVar(name="WANDB_RUN_GROUP", value=self.wandb_run_group),
        ]
        spec = beaker.ExperimentSpec(
            budget="ai2/climate",
            description="Do inference with ACE-SOM.",
            tasks=[
                beaker.TaskSpec(
                    name=self.name,
                    image=beaker.ImageSource(beaker=self.image),
                    command=[
                        "python",
                        "-m",
                        self.entrypoint,
                        f"{DATASET_CONFIG_MOUNTPATH}/{DATASET_CONFIG_FILENAME}",
                    ],
                    result=beaker.ResultSpec(path="/output"),
                    resources=beaker.TaskResources(gpu_count=1, shared_memory="50GiB"),
                    context=beaker.TaskContext(priority=self.priority, preemptible=True),
                    constraints=beaker.Constraints(cluster=CLUSTERS),
                    env_vars=env_vars,
                    datasets=self.datasets,
                )
            ],
        )
        return spec

    def submit(self):
        spec = self.experiment_spec
        try:
            experiment = self.client.experiment.create(self.name, spec)
            print(
                f"Experiment {self.name} created. See https://beaker.org/ex/{experiment.id}"
            )
        except beaker.exceptions.ExperimentConflict:
            print(
                f"Failed to create experiment {self.name} because it already exists. "
                "Skipping experiment creation. If you want to submit this experiment, "
                "delete the existing experiment with the same name, or rename the new "
                "experiment."
            )


@dataclasses.dataclass
class EvaluatorRollout(BaseRollout):
    client: beaker.Beaker
    image: str
    model: Model
    steps: int
    reference_dataset: ReferenceDataset
    prediction_dataset: Optional[ReferenceDataset] = None
    tag: Optional[str] = None
    overlay_template: Optional[Dict] = None
    wandb_run_group: str = "ace2-som"
    priority: str = "normal"

    @property
    def name(self):
        date = str(datetime.date.today())
        if self.tag is None:
            prefix = date
        else:
            prefix = "-".join([date, self.tag])
        if self.prediction_dataset is None:
            return "-".join([prefix, self.reference_dataset.name, self.model.name])
        else:
            return "-".join(
                [
                    prefix,
                    "data-only",
                    self.prediction_dataset.name,
                    self.model.name,
                ]
            )

    @property
    def overlay(self):
        if self.overlay_template is None:
            result = copy.deepcopy(DEFAULT_OVERLAY_TEMPLATE)
        else:
            result = copy.deepcopy(self.overlay_template)
        result["n_forward_steps"] = self.steps
        result["loader"]["dataset"]["data_path"] = self.reference_dataset.path
        result["loader"]["start_indices"] = {
            "times": self.reference_dataset.start_times
        }
        if self.prediction_dataset is not None:
            result["prediction_loader"] = {
                "dataset": {"data_path": self.prediction_dataset.path},
                "start_indices": {"times": self.prediction_dataset.start_times},
                "num_data_workers": 8,
            }
        return result

    @property
    def datasets(self):
        return self.standard_datasets

    @property
    def entrypoint(self):
        return EVALUATOR


@dataclasses.dataclass
class InitialCondition:
    name: str
    path: str
    time: str
    data_source: Optional[beaker.DataSource] = None


@dataclasses.dataclass
class ForcingDataset:
    name: str
    path: str


@dataclasses.dataclass
class InferenceRollout(BaseRollout):
    client: beaker.Beaker
    image: str
    model: Model
    steps: int
    forcing_dataset: ForcingDataset
    initial_condition: InitialCondition
    tag: str
    overlay_template: Optional[Dict] = None
    wandb_run_group: str = "ace2-som"
    priority: str = "normal"

    @property
    def name(self):
        date = str(datetime.date.today())
        return "-".join(
            [
                date,
                self.tag,
                self.initial_condition.name,
                self.forcing_dataset.name,
                self.model.name,
            ]
        )

    @property
    def overlay(self):
        if self.overlay_template is None:
            result = copy.deepcopy(DEFAULT_OVERLAY_TEMPLATE)
        else:
            result = copy.deepcopy(self.overlay_template)
        result["n_forward_steps"] = self.steps
        result["forcing_loader"]["dataset"]["data_path"] = self.forcing_dataset.path
        result["initial_condition"]["path"] = self.initial_condition.path
        result["initial_condition"]["start_indices"] = {
            "times": [self.initial_condition.time]
        }
        return result

    @property
    def datasets(self):
        result = self.standard_datasets
        if self.initial_condition.data_source is not None:
            result.append(
                beaker.DataMount(
                    mount_path=INITIAL_CONDITION_MOUNT_PATH,
                    source=self.initial_condition.data_source,
                )
            )
        return result

    @property
    def entrypoint(self):
        return INFERENCE
