rule initialize_zarr_regridded_tape:
    input:
        kerchunk_json=rules.kerchunk_regridded_tape.output.json
    output:
        store=directory(OUTPUT_DIRECTORY / "{target_grid}" / "{run}" / "{tape}.zarr")
    resources:
        mem="64G",
        partition="dtn_f5_f6",
        constraint="f5"
    script:
        "../scripts/initialize_store.py"


rule write_partition_zarr_regridded_tape:
    input:
        kerchunk_json=rules.kerchunk_regridded_tape.output.json,
        store=rules.initialize_zarr_regridded_tape.output.store
    output:
        touch(f"xpartition-sentinels/{CONFIG}/{{target_grid}}/{{run}}/{{tape}}.{{rank}}.out")
    resources:
        mem="64G",
        partition="dtn_f5_f6",
        constraint="f5"
    params:
        ranks=config["xpartition_ranks"],
        rank="{rank}",
    wildcard_constraints:
        rank="\d+"
    script:
        "../scripts/write_partition.py"


rule zarr_regridded_tape:
    input:
        kerchunk_json=rules.kerchunk_regridded_tape.output.json,
        sentinels=expand(
            rules.write_partition_zarr_regridded_tape.output,
            rank=XPARTITION_RANKS,
            run="{run}",
            tape="{tape}",
            target_grid="{target_grid}"
        )
    output:
        touch(f"xpartition-sentinels/{CONFIG}/{{target_grid}}/{{run}}/{{tape}}.complete.out")


rule zarr_scalar_tape:
    input:
        kerchunk_json=rules.kerchunk_scalar_tape.output.json,
    output:
        store=directory(OUTPUT_DIRECTORY / "{target_grid}" / "{run}" / "{scalar_tape}.zarr")
    resources:
        mem="64G",
        partition="dtn_f5_f6",
        constraint="f5"
    script:
        "../scripts/zarr_scalar_tape.py"
