import argparse

import cftime
import fv3config
import xarray as xr
import yaml

from typing import Dict, Tuple


UNITS = {"seconds", "minutes", "hours", "days", "months"}
NON_MONTH_UNITS = {"seconds", "minutes", "hours", "days"}
DATE_FORMAT = "%Y%m%d%H"


def get_duration_config(config: dict) -> Dict[str, int]:
    coupler_nml = config["namelist"].get("coupler_nml", {})
    return {unit: coupler_nml.get(unit, 0) for unit in UNITS}


def duration_timedelta_representable(config: dict) -> bool:
    duration_config = get_duration_config(config)
    return duration_config["months"] == 0


def duration_months_representable(config: dict) -> bool:
    duration_config = get_duration_config(config)
    return all(duration_config[unit] == 0 for unit in NON_MONTH_UNITS)


def get_initial_date(config: dict) -> cftime.DatetimeJulian:
    initial_date_tuple = config["namelist"]["coupler_nml"]["current_date"]
    return cftime.DatetimeJulian(*initial_date_tuple)


def get_date(config: dict, segment: int, end_date: bool = False):
    """Return the date associated with the start or end of the segment

    Args:
        config: dict
            fv3config configuration dictionary for the simulation
        segment: int
            segment of the simulation (note this is zero indexed; the first
            segment is segment zero).
        end_date: bool
            whether to return the end date of the segment instead of the start

    Returns:
        date: cftime.DatetimeJulian
    """
    initial_date = get_initial_date(config)

    if duration_timedelta_representable(config):
        duration = fv3config.get_run_duration(config)
        if end_date:
            segment = segment + 1
        date = initial_date + segment * duration
    elif duration_months_representable(config) and initial_date.day == 1:
        duration = get_duration_config(config)["months"]
        if end_date:
            segment = segment + 1
        freq = f"{duration}MS"
        periods = segment + 1  # segment == 1 to return initial date
        date_range = xr.cftime_range(initial_date, freq=freq, periods=periods)
        date = date_range.values[-1]
    else:
        raise ValueError("Segment length and initial date combination not supported.")

    return date


def get_segment_dates(config: dict, segment_number: int) -> Tuple[str, str]:
    start_date = get_date(config, segment_number, end_date=False)
    end_date = get_date(config, segment_number, end_date=True)
    return start_date.strftime(DATE_FORMAT), end_date.strftime(DATE_FORMAT)
