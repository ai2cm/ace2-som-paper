rule archive_diagnostics_file:
    input:
        get_lustre_history_file_path
    output:
        touch(f"archive-sentinels/{CONFIG}/archived-diagnostics-{{run}}-{{segment}}-{{tape}}-{{tile}}.out")
    resources:
        mem="8G",
        partition="dtn_f5_f6"
    params:
        archive_directory=get_archive_history_file_directory
    shell:
        """
        module load gcp
        gcp -cd {input} gfdl:{params.archive_directory}/
        """


rule archive_scalar_diagnostics_file:
    input:
        get_lustre_scalar_history_file_path
    output:
        touch(f"archive-sentinels/{CONFIG}/archived-scalar-diagnostics-{{run}}-{{segment}}-{{scalar_tape}}.out")
    resources:
        mem="8G",
        partition="dtn_f5_f6"
    params:
        archive_directory=get_archive_history_file_directory
    shell:
        """
        module load gcp
        gcp -cd {input} gfdl:{params.archive_directory}/
        """


rule tar_restart_directory:
    input:
        source=get_lustre_restart_directory
    output:
        destination=WORKING_DIRECTORY / "{run}" / "restart" / "{restart_segment}.tar"
    resources:
        mem="8G",
        partition="dtn_f5_f6",
        constraint="f5"
    shell:
        """tar -cvf {output.destination} -C {input.source} ."""


rule archive_restart_directory:
    input:
        rules.tar_restart_directory.output.destination
    output:
        touch(f"archive-sentinels/{CONFIG}/archived-restarts-{{run}}-{{restart_segment}}.out")
    resources:
        mem="8G",
        partition="dtn_f5_f6",
        constraint="f5"
    params:
        archive_directory=get_archive_restart_file_directory
    shell:
        """
        module load gcp
        gcp -cd {input} gfdl:{params.archive_directory}/
        """


rule tar_ascii_directory:
    input:
        source=get_lustre_ascii_directory
    output:
        destination=WORKING_DIRECTORY / "{run}" / "ascii" / "{segment}.tar"
    resources:
        mem="8G",
        partition="dtn_f5_f6",
        constraint="f5"
    shell:
        """tar -cvf {output.destination} -C {input.source} ."""


rule archive_ascii_directory:
    input:
        rules.tar_ascii_directory.output.destination
    output:
        touch(f"archive-sentinels/{CONFIG}/archived-ascii-{{run}}-{{segment}}.out")
    resources:
        mem="8G",
        partition="dtn_f5_f6",
        constraint="f5"
    params:
        archive_directory=get_archive_ascii_file_directory
    shell:
        """
        module load gcp
        gcp -cd {input} gfdl:{params.archive_directory}/
        """


rule archive_diagnostics_file_all:
    input:
        [
            expand(
                rules.archive_diagnostics_file.output,
                run=run,
                segment=SEGMENTS[run],
                tape=TAPES[run],
                tile=TILES,
            )
        for run in RUNS],
        [
            expand(
                rules.archive_scalar_diagnostics_file.output,
                run=run,
                segment=SEGMENTS[run],
                scalar_tape=SCALAR_TAPES[run],
            )
        for run in RUNS]


rule archive_restart_directory_all:
    input:
        [
            expand(
                rules.archive_restart_directory.output,
                run=run,
                restart_segment=RESTART_SEGMENTS[run],
                tape=TAPES[run]
            )
        for run in RUNS]


rule archive_ascii_directory_all:
    input:
        [
            expand(
                rules.archive_ascii_directory.output,
                run=run,
                segment=SEGMENTS[run],
                tape=TAPES[run]
            )
        for run in RUNS]
        
