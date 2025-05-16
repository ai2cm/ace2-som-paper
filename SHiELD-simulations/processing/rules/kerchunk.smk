rule kerchunk_regridded_segment_tape:
    input:
        file=rules.regrid_tape.output.file
    output:
        json=WORKING_DIRECTORY / "{target_grid}" / "{run}" / "{segment}" / "kerchunk-single-file-{tape}.json"
    resources:
        mem="8G",
        partition="dtn_f5_f6",
        constraint="f5"
    script:
        "../scripts/kerchunk_single_file.py"


rule kerchunk_regridded_tape:
    input:
        kerchunk_regridded_tape_input
    output:
        json=WORKING_DIRECTORY / "{target_grid}" / "{run}" / "kerchunk-combined-{tape}.json"
    resources:
        mem="16G",
        partition="dtn_f5_f6",
        constraint="f5"
    script:
        "../scripts/kerchunk_multi_file.py"


rule kerchunk_scalar_segment_tape:
    input:
        file=kerchunk_scalar_segment_tape_input
    output:
        json=WORKING_DIRECTORY / "{target_grid}" / "{run}" / "{segment}" / "kerchunk-single-scalar-file-{scalar_tape}.json"
    resources:
        mem="8G",
        partition="dtn_f5_f6",
        constraint="f5"
    script:
        "../scripts/kerchunk_single_file.py"


rule kerchunk_scalar_tape:
    input:
        file=kerchunk_scalar_tape_input
    output:
        json=WORKING_DIRECTORY / "{target_grid}" / "{run}" / "kerchunk-scalar-combined-{scalar_tape}.json"
    resources:
        mem="8G",
        partition="dtn_f5_f6",
        constraint="f5"
    script:
        "../scripts/kerchunk_multi_file.py"
