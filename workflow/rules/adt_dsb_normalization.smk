import os

WORKDIR = os.getcwd()


# Rule normalizing the ADT counts using R package "dsb"
# Input is cellranger ADT output
# Output is SCE object in RDS file
rule dsb_normalize_adt:
    input:
        CellrangerADT="results/cellranger_adt/{sample}.matrix.mtx",
    output:
        normalized="results/dsb_normalize_adt/{sample}.dsb_normalize_adt.RDS",
    singularity:
        config["tools"]["dsb_normalize_adt"]["singularity"]
    params:
        adt_raw_dir=WORKDIR
        + "/results/cellranger_adt/{sample}/outs/raw_feature_bc_matrix/",
        adt_filtered_dir=WORKDIR
        + "/results/cellranger_adt/{sample}/outs/filtered_feature_bc_matrix/",
        outdir="results/dsb_normalize_adt/",
    threads: config["computingResources"]["threads"]["high"]
    log:
        "logs/dsb_normalize_adt/{sample}.log",
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["high"],
        runtime=config["computingResources"]["runtime"]["high"],
    benchmark:
        "results/dsb_normalize_adt/benchmark/{sample}.benchmark"
    shell:
        "Rscript workflow/scripts/dsb_normalize_adt.R "
        "--adt_raw_dir {params.adt_raw_dir} "
        "--adt_filtered_dir {params.adt_filtered_dir} "
        "--outdir {params.outdir} "
        "--sample {wildcards.sample} "
