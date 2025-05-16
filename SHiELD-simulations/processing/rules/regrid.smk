rule remap_file:
    output:
        remap_file=WORKING_DIRECTORY / "remap_{remap_source_grid}_to_{remap_target_grid}.nc"
    resources:
        mem="32G",
        partition="dtn_f5_f6",
        constraint="f5"
    envmodules:
        "fre-nctools/2023.01"
    script:
        "../scripts/generate_remap_file.py"


# Note we ensure we generate ALL remap files ahead of time, since the source
# grid resolution is determined dynamically. This only needs to be done once
# anyway, so it is good to get out of the way.
rule regrid_tape:
    input:
        expand(rules.remap_file.output, remap_source_grid=SOURCE_GRIDS, remap_target_grid=TARGET_GRIDS),
        regrid_tape_input
    output:
        file=WORKING_DIRECTORY / "{target_grid}" / "{run}" / "{segment}" / "{tape}.nc"
    resources:
        mem="32G",
        partition="dtn_f5_f6",
        constraint="f5"
    params:
        segment_directory=get_segment_directory
    envmodules:
        "fre-nctools/2023.01"
    script:
        "../scripts/regrid.py"
