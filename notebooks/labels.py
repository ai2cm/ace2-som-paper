import unicodedata

# mapping from text name to symbol
air_temp = {f'air_temperature_{k}': f'$T_{k}$' for k in range(8)}
east_wind = {f'eastward_wind_{k}': f'$u_{k}$' for k in range(8)}
north_wind = {f'northward_wind_{k}': f'$v_{k}$' for k in range(8)}
water = {f'specific_total_water_{k}': f'$q^T_{k}$' for k in range(8)}
prognostic_1d = {
    'surface_temperature': '$T_s$',
    'PRESsfc': '$p_s$',
    'total_water_path': 'TWP'
}
prognostic = {**prognostic_1d, **air_temp, **water, **east_wind, **north_wind}
forcing = {'HGTsfc': '$z_s$', 'DSWRFtoa': 'DSWRFtoa'}
diagnostic = {
    'PRATEsfc': '$P$',
    'LHTFLsfc': '$LHF$',
    'SHTFLsfc': '$SHF$',
    'tendency_of_total_water_path_due_to_advection': 'TWP adv tend',
    'USWRFtoa': 'USWRFtoa',
    'ULWRFtoa': 'ULWRFtoa',
    'DSWRFsfc': 'DSWRFsfc',
    'DLWRFsfc': 'DLWRFsfc',
    'USWRFsfc': 'USWRFsfc',
    'ULWRFsfc': 'ULWRFsfc',
}
PROGNOSTIC_VARIABLES = list(prognostic)
DIAGNOSTIC_VARIABLES = list(diagnostic)
LONG_NAMES = {**prognostic, **forcing, **diagnostic}

# mapping from text name to symbol using unicode only (necessary for altair!)
digits = ['Zero', 'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven']
air_temp = {f'air_temperature_{k}': 'T'+unicodedata.lookup(f'Subscript {d}') for k, d in enumerate(digits)}
east_wind = {f'eastward_wind_{k}': 'u'+unicodedata.lookup(f'Subscript {d}') for k, d in enumerate(digits)}
north_wind = {f'northward_wind_{k}': 'v'+unicodedata.lookup(f'Subscript {d}') for k, d in enumerate(digits)}
water = {f'specific_total_water_{k}': 'q\u1D40'+unicodedata.lookup(f'Subscript {d}') for k, d in enumerate(digits)}
prognostic_1d = {
    'surface_temperature': 'T\u209b',
    'PRESsfc': 'p\u209b',
    'total_water_path': 'TWP',
}
prognostic = {**prognostic_1d, **air_temp, **water, **east_wind, **north_wind}
forcing = {'HGTsfc': 'z\u209b', 'DSWRFtoa': 'DSWRFtoa'}
diagnostic = {
    'PRATEsfc': 'P',
    'LHTFLsfc': 'LHF',
    'SHTFLsfc': 'SHF',
    'tendency_of_total_water_path_due_to_advection': 'TWP adv tend',
    'USWRFtoa': 'USWRFtoa',
    'ULWRFtoa': 'ULWRFtoa',
    'DSWRFsfc': 'DSWRFsfc',
    'DLWRFsfc': 'DLWRFsfc',
    'USWRFsfc': 'USWRFsfc',
    'ULWRFsfc': 'ULWRFsfc',
    'Q2m': 'Q2m',
    'TMP2m': 'T2m',
    'UGRD10m': 'u10m',
    'VGRD10m': 'v10m',
    'TMP850': "T850",
    "h500": 'h500'
}
UNICODE_NAMES = {**prognostic, **forcing, **diagnostic}

# mapping from text name to units (latex and unicode-only)
air_temp = {f'air_temperature_{k}': 'K' for k in range(7, -1, -1)}
east_wind = {f'eastward_wind_{k}': 'm/s' for k in range(7, -1, -1)}
north_wind = {f'northward_wind_{k}': 'm/s' for k in range(7, -1, -1)}
water = {f'specific_total_water_{k}': 'kg/kg' for k in range(7, -1, -1)}
prognostic_1d = {
    'surface_temperature': 'K',
    'PRESsfc': 'Pa',
    'total_water_path': 'mm',
}
forcing = {'HGTsfc': 'm', 'DSWRFtoa': 'W/m$^2$'}
diagnostic = {
    'PRATEsfc': 'mm/day',
    'LHTFLsfc': 'W/m$^2$',
    'SHTFLsfc': 'W/m$^2$',
    'tendency_of_total_water_path_due_to_advection': 'kg/m$^2$/s',
    'USWRFtoa': 'W/m$^2$',
    'ULWRFtoa': 'W/m$^2$',
    'DSWRFsfc': 'W/m$^2$',
    'DLWRFsfc': 'W/m$^2$',
    'USWRFsfc': 'W/m$^2$',
    'ULWRFsfc': 'W/m$^2$',
    'Q2m': 'kg/kg',
    'TMP2m': 'K',
    'UGRD10m': 'm/s',
    'VGRD10m': 'm/s',
    'TMP850': "K",
    "h500": 'm'
}
diagnostic_unicode = {
    'PRATEsfc': 'kg/m\u00b2/s',
    'LHTFLsfc': 'W/m\u00b2',
    'SHTFLsfc': 'W/m\u00b2',
    'tendency_of_total_water_path_due_to_advection': 'kg/m\u00b2/s',
    'USWRFtoa': 'W/m\u00b2',
    'ULWRFtoa': 'W/m\u00b2',
    'DSWRFsfc': 'W/m\u00b2',
    'DLWRFsfc': 'W/m\u00b2',
    'USWRFsfc': 'W/m\u00b2',
    'ULWRFsfc': 'W/m\u00b2',
    "column_moist_static_energy": "J/m\u00b2",
    'Q2m': 'kg/kg',
    'TMP2m': 'K',
    'UGRD10m': 'm/s',
    'VGRD10m': 'm/s',
    'TMP850': "K",
    "h500": 'm'
}
UNITS = {**prognostic_1d, **air_temp, **water, **east_wind, **north_wind, **forcing, **diagnostic}
UNITS_UNICODE = {**prognostic_1d, **air_temp, **water, **east_wind, **north_wind, **forcing, **diagnostic_unicode}
