import cftime
import pytest


from shield_run.dates import (
    duration_months_representable,
    duration_timedelta_representable,
    get_date,
    get_duration_config,
    get_segment_dates,
)


def test_get_duration_config():
    config = {"namelist": {"coupler_nml": {"days": 1, "seconds": 5}}}
    expected = {"months": 0, "days": 1, "hours": 0, "minutes": 0, "seconds": 5}
    result = get_duration_config(config)
    assert result == expected


@pytest.mark.parametrize(("months", "expected"), [(0, True), (1, False)])
def test_duration_timedelta_representable(months, expected):
    config = {"namelist": {"coupler_nml": {"months": months}}}
    assert duration_timedelta_representable(config) == expected


@pytest.mark.parametrize("units", ["seconds", "minutes", "hours", "days"])
def test_duration_months_representable(units):
    config = {"namelist": {"coupler_nml": {units: 1}}}
    assert duration_months_representable(config) == False


def test_get_date_timedelta_representable():
    config = {
        "namelist": {
            "coupler_nml": {
                "current_date": [2020, 2, 1, 0, 0, 0],
                "days": 1,
                "hours": 10,
            }
        }
    }
    segment = 1

    result_start = get_date(config, segment, end_date=False)
    expected_start = cftime.DatetimeJulian(2020, 2, 2, 10)
    assert result_start == expected_start

    result_end = get_date(config, segment, end_date=True)
    expected_end = cftime.DatetimeJulian(2020, 2, 3, 20)
    assert result_end == expected_end


def test_get_date_months_representable():
    config = {
        "namelist": {
            "coupler_nml": {"current_date": [2020, 2, 1, 0, 0, 0], "months": 2}
        }
    }
    segment = 1

    result_start = get_date(config, segment, end_date=False)
    expected_start = cftime.DatetimeJulian(2020, 4, 1)
    assert result_start == expected_start

    result_end = get_date(config, segment, end_date=True)
    expected_end = cftime.DatetimeJulian(2020, 6, 1)
    assert result_end == expected_end


def test_get_date_invalid_months_representable():
    config = {
        "namelist": {
            "coupler_nml": {"current_date": [2020, 2, 5, 0, 0, 0], "months": 2}
        }
    }
    segment = 1
    with pytest.raises(ValueError, match="Segment length"):
        get_date(config, segment)


def test_get_segment_dates():
    config = {
        "namelist": {
            "coupler_nml": {
                "current_date": [2020, 2, 1, 0, 0, 0],
                "days": 1,
                "hours": 10,
            }
        }
    }
    segment = 1
    expected = ("2020020210", "2020020320")
    result = get_segment_dates(config, segment)
    assert result == expected
