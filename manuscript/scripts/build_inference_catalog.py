import logging

import beaker
import pandas as pd
import wandb


from typing import Dict, Tuple


CLIMATES = ["1xCO2", "2xCO2", "3xCO2", "4xCO2"]
INITIAL_CONDITIONS = [f"ic{n}" for n in range(1, 6)]
SEEDS = 4


def beaker_dataset_id_from_name(name):
    client = beaker.Beaker.from_env()
    experiment_client = beaker.services.ExperimentClient(client)
    return (
        experiment_client
        .get(name)
        .jobs[-1]
        .result
        .beaker
    )


def get_shield_experiment_name(model, climate, initial_condition):
    return f"2024-11-24-data-only-{initial_condition}-{climate}-{model}"


def get_ace_experiment_name(model, climate, initial_condition):
    return f"2024-11-24-ensemble-member-{initial_condition}-{climate}-{model}"


def get_shield_experiment_tuned_cdmbgwd_name_extreme_PRATEsfc(
    model, climate, initial_condition
):
    return (
        f"2024-11-26-extreme-PRATEsfc-data-only-{initial_condition}-{climate}-{model}"
    )


def get_shield_experiment_name_extreme_PRATEsfc(model, climate, initial_condition):
    return (
        f"2024-11-25-extreme-PRATEsfc-data-only-{initial_condition}-{climate}-{model}"
    )


def get_ace_experiment_name_extreme_PRATEsfc(model, climate, initial_condition):
    return f"2024-11-25-ensemble-member-extreme-PRATEsfc-{initial_condition}-{climate}-{model}"


GET_EQUILIBRIUM_NAME_MAPPING = {
    "SHiELD-SOM-C24-tuned-cdmbgwd": get_shield_experiment_name,
    "SHiELD-SOM-C96": get_shield_experiment_name,
    **{f"ACE2-SOM-multi-climate-RS{s}": get_ace_experiment_name for s in range(SEEDS)},
    **{f"ACE2-SOM-increasing-co2-RS{s}": get_ace_experiment_name for s in range(SEEDS)},
}


GET_INCREASING_CO2_NAME_MAPPING = {
    "SHiELD-SOM-C24-tuned-cdmbgwd": "2024-11-24-data-only-increasing-CO2-SHiELD-SOM-C24-tuned-cdmbgwd",
    "SHiELD-SOM-C96": "2024-11-24-data-only-increasing-CO2-SHiELD-SOM-C96",
    **{
        f"ACE2-SOM-multi-climate-RS{s}": f"2024-11-24-increasing-co2-ACE2-SOM-multi-climate-RS{s}"
        for s in range(SEEDS)
    },
    **{
        f"ACE2-SOM-increasing-co2-RS{s}": f"2024-11-24-increasing-co2-ACE2-SOM-increasing-co2-RS{s}"
        for s in range(SEEDS)
    },
}


GET_EQUILIBRIUM_NAME_MAPPING_EXTREME_PRATEsfc = {
    "SHiELD-SOM-C24-tuned-cdmbgwd": get_shield_experiment_tuned_cdmbgwd_name_extreme_PRATEsfc,
    "SHiELD-SOM-C96": get_shield_experiment_name_extreme_PRATEsfc,
    **{
        f"ACE2-SOM-multi-climate-RS{s}": get_ace_experiment_name_extreme_PRATEsfc
        for s in range(SEEDS)
    },
    **{
        f"ACE2-SOM-increasing-co2-RS{s}": get_ace_experiment_name_extreme_PRATEsfc
        for s in range(SEEDS)
    },
}

GET_ABRUPT_4xCO2_MAPPING = {
    "SHiELD-SOM-C24-tuned-cdmbgwd": "2024-11-24-data-only-abrupt-4xCO2-SHiELD-SOM-C24-tuned-cdmbgwd",
    "SHiELD-SOM-C96": "2024-11-24-data-only-abrupt-4xCO2-SHiELD-SOM-C96",
    **{
        f"ACE2-SOM-multi-climate-RS{s}": f"2024-11-24-abrupt-4xCO2-ACE2-SOM-multi-climate-RS{s}"
        for s in range(SEEDS)
    },
    **{
        f"ACE2-SOM-increasing-co2-RS{s}": f"2024-11-24-abrupt-4xCO2-ACE2-SOM-increasing-CO2-RS{s}"
        for s in range(SEEDS)
    },
}

GET_ABRUPT_4xCO2_ABBREVIATED_MAPPING = {
    "SHiELD-SOM-C24-tuned-cdmbgwd": "2024-11-24-abbreviated-data-only-abrupt-4xCO2-SHiELD-SOM-C24-tuned-cdmbgwd",
    "SHiELD-SOM-C96": "2024-11-24-abbreviated-data-only-abrupt-4xCO2-SHiELD-SOM-C96",
    **{
        f"ACE2-SOM-multi-climate-RS{s}": f"2024-11-24-abbreviated-abrupt-4xCO2-ACE2-SOM-multi-climate-RS{s}"
        for s in range(SEEDS)
    },
    **{
        f"ACE2-SOM-increasing-co2-RS{s}": f"2024-11-24-abbreviated-abrupt-4xCO2-ACE2-SOM-increasing-CO2-RS{s}"
        for s in range(SEEDS)
    },
}

GET_RADIATION_MULTI_CALLS_MAPPING = {
    **{
        f"ACE2-SOM-multi-climate-RS{s}": f"2024-11-30-radiation-multi-calls-ACE2-SOM-multi-climate-RS{s}"
        for s in range(SEEDS)
    },
    **{
        f"ACE2-SOM-increasing-co2-RS{s}": f"2024-11-30-radiation-multi-calls-ACE2-SOM-increasing-CO2-RS{s}"
        for s in range(SEEDS)
    },
}


def wandb_id_from_name(name):
    api = wandb.Api()
    runs = api.runs("ai2cm/ace-som", filters={"display_name": {"$regex": name}})
    completed = [run for run in runs if run.state == "finished"]
    if len(completed) > 1:
        raise ValueError(f"Ambiguous name {name}")
    elif len(completed) == 0:
        logging.info(f"No WandB run found for {name}")
        return None
    else:
        (result,) = completed
        return "/".join(result.path)


def fetch_equilibrium_ids(model: str) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
    dataset_ids = {}
    get_name = GET_EQUILIBRIUM_NAME_MAPPING[model]
    for climate in CLIMATES:
        for initial_condition in INITIAL_CONDITIONS:
            key = model, climate, initial_condition, None
            name = get_name(model, climate, initial_condition)
            beaker_id = beaker_dataset_id_from_name(name)
            wandb_id = wandb_id_from_name(name)
            ids = beaker_id, wandb_id
            dataset_ids[key] = ids
            logging.info(f"Fetched {key}: {ids}")
    return dataset_ids


def fetch_increasing_co2_ids(model: str) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
    key = model, "increasing-CO2", None, None
    name = GET_INCREASING_CO2_NAME_MAPPING[model]
    beaker_id = beaker_dataset_id_from_name(name)
    wandb_id = wandb_id_from_name(name)
    ids = beaker_id, wandb_id
    logging.info(f"Fetched {key}: {ids}")
    return {key: ids}


def fetch_equilibrium_ids_extreme_precipitation(
    model: str,
) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
    dataset_ids = {}
    get_name = GET_EQUILIBRIUM_NAME_MAPPING_EXTREME_PRATEsfc[model]
    for climate in CLIMATES:
        for initial_condition in INITIAL_CONDITIONS:
            key = model, climate, initial_condition, "extreme-precipitation"
            name = get_name(model, climate, initial_condition)
            beaker_id = beaker_dataset_id_from_name(name)
            wandb_id = wandb_id_from_name(name)
            ids = beaker_id, wandb_id
            dataset_ids[key] = ids
            logging.info(f"Fetched {key}: {ids}")
    return dataset_ids


def fetch_abrupt_4xCO2_ids(model: str) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
    key = model, "abrupt-4xCO2", None, None
    name = GET_ABRUPT_4xCO2_MAPPING[model]
    beaker_id = beaker_dataset_id_from_name(name)
    wandb_id = wandb_id_from_name(name)
    ids = beaker_id, wandb_id
    logging.info(f"Fetched {key}: {ids}")
    return {key: ids}


def fetch_abrupt_4xCO2_abbreviated_ids(
    model: str,
) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
    key = model, "abrupt-4xCO2", None, "abbreviated"
    name = GET_ABRUPT_4xCO2_ABBREVIATED_MAPPING[model]
    beaker_id = beaker_dataset_id_from_name(name)
    wandb_id = wandb_id_from_name(name)
    ids = beaker_id, wandb_id
    logging.info(f"Fetched {key}: {ids}")
    return {key: ids}


def fetch_radiation_multi_calls_ids(
    model: str,
) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
    key = model, "radiation-multi-calls", None, None
    name = GET_RADIATION_MULTI_CALLS_MAPPING[model]
    beaker_id = beaker_dataset_id_from_name(name)
    wandb_id = wandb_id_from_name(name)
    ids = beaker_id, wandb_id
    logging.info(f"Fetched {key}: {ids}")
    return {key: ids}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    equilibrium_climate_inference = {
        **fetch_equilibrium_ids("SHiELD-SOM-C24-tuned-cdmbgwd"),
        **fetch_equilibrium_ids("SHiELD-SOM-C96"),
        **fetch_equilibrium_ids("ACE2-SOM-multi-climate-RS0"),
        **fetch_equilibrium_ids("ACE2-SOM-multi-climate-RS1"),
        **fetch_equilibrium_ids("ACE2-SOM-multi-climate-RS2"),
        **fetch_equilibrium_ids("ACE2-SOM-multi-climate-RS3"),
        **fetch_equilibrium_ids("ACE2-SOM-increasing-co2-RS0"),
        **fetch_equilibrium_ids("ACE2-SOM-increasing-co2-RS1"),
        **fetch_equilibrium_ids("ACE2-SOM-increasing-co2-RS2"),
        **fetch_equilibrium_ids("ACE2-SOM-increasing-co2-RS3"),
    }
    increasing_co2_inference = {
        **fetch_increasing_co2_ids("SHiELD-SOM-C24-tuned-cdmbgwd"),
        **fetch_increasing_co2_ids("SHiELD-SOM-C96"),
        **fetch_increasing_co2_ids("ACE2-SOM-multi-climate-RS0"),
        **fetch_increasing_co2_ids("ACE2-SOM-multi-climate-RS1"),
        **fetch_increasing_co2_ids("ACE2-SOM-multi-climate-RS2"),
        **fetch_increasing_co2_ids("ACE2-SOM-multi-climate-RS3"),
        **fetch_increasing_co2_ids("ACE2-SOM-increasing-co2-RS0"),
        **fetch_increasing_co2_ids("ACE2-SOM-increasing-co2-RS1"),
        **fetch_increasing_co2_ids("ACE2-SOM-increasing-co2-RS2"),
        **fetch_increasing_co2_ids("ACE2-SOM-increasing-co2-RS3"),
    }
    equilibrium_climate_inference_extreme_precipitation = {
        **fetch_equilibrium_ids_extreme_precipitation("SHiELD-SOM-C24-tuned-cdmbgwd"),
        **fetch_equilibrium_ids_extreme_precipitation("SHiELD-SOM-C96"),
        **fetch_equilibrium_ids_extreme_precipitation("ACE2-SOM-multi-climate-RS3"),
    }
    abrupt_4xCO2_inference = {
        **fetch_abrupt_4xCO2_ids("SHiELD-SOM-C24-tuned-cdmbgwd"),
        **fetch_abrupt_4xCO2_ids("SHiELD-SOM-C96"),
        **fetch_abrupt_4xCO2_ids("ACE2-SOM-multi-climate-RS3"),
    }
    abrupt_4xCO2_inference_abbreviated = {
        **fetch_abrupt_4xCO2_abbreviated_ids("SHiELD-SOM-C24-tuned-cdmbgwd"),
        **fetch_abrupt_4xCO2_abbreviated_ids("SHiELD-SOM-C96"),
        **fetch_abrupt_4xCO2_abbreviated_ids("ACE2-SOM-multi-climate-RS3"),
    }
    radiation_multi_calls = {
        **fetch_radiation_multi_calls_ids("ACE2-SOM-multi-climate-RS0"),
        **fetch_radiation_multi_calls_ids("ACE2-SOM-multi-climate-RS1"),
        **fetch_radiation_multi_calls_ids("ACE2-SOM-multi-climate-RS2"),
        **fetch_radiation_multi_calls_ids("ACE2-SOM-multi-climate-RS3"),
        **fetch_radiation_multi_calls_ids("ACE2-SOM-increasing-co2-RS0"),
        **fetch_radiation_multi_calls_ids("ACE2-SOM-increasing-co2-RS1"),
        **fetch_radiation_multi_calls_ids("ACE2-SOM-increasing-co2-RS2"),
        **fetch_radiation_multi_calls_ids("ACE2-SOM-increasing-co2-RS3"),
    }
    inference_runs = {
        **equilibrium_climate_inference,
        **increasing_co2_inference,
        **equilibrium_climate_inference_extreme_precipitation,
        **abrupt_4xCO2_inference,
        **abrupt_4xCO2_inference_abbreviated,
        **radiation_multi_calls,
    }
    index = pd.MultiIndex.from_tuples(
        inference_runs.keys(), names=["model", "forcing", "initial_condition", "tag"]
    )
    df = pd.DataFrame(
        inference_runs.values(), columns=["beaker_id", "wandb_id"], index=index
    )
    df = df.reset_index()
    print(df)
    df.to_csv("inference-run-catalog.csv", index=False)
