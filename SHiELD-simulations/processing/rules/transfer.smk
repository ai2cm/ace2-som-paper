rule transfer_tape:
    input:
        sentinel=rules.zarr_regridded_tape.output,
        store=rules.initialize_zarr_regridded_tape.output.store
    output:
        touch(f"transfer-sentinels/{CONFIG}/regridded-zarrs/{{target_grid}}/{{run}}/{{tape}}.out")
    resources:
        partition="dtn_f5_f6",
        mem="8GB",
        constraint="f5"
    params:
        gcs_path=f"{GCS_BUCKET}/regridded-zarrs/{{target_grid}}/{{run}}/"
    retries: 5
    shell:
        """
        gsutil -m cp -n -r {input.store} {params.gcs_path}
        """


rule transfer_scalar_tape:
    input:
        store=rules.zarr_scalar_tape.output.store
    output:
        touch(f"transfer-sentinels-scalar/{CONFIG}/regridded-zarrs/{{target_grid}}/{{run}}/{{scalar_tape}}.out")
    resources:
        partition="dtn_f5_f6",
        mem="8GB",
        constraint="f5"
    params:
        gcs_path=f"{GCS_BUCKET}/regridded-zarrs/{{target_grid}}/{{run}}/"
    retries: 5
    shell:
        """
        gsutil -m cp -n -r {input.store} {params.gcs_path}
        """
