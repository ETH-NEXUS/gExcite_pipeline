import os

WORKDIR = os.getcwd()


# Cellranger call to process the raw GEX samples
# input: fastq directory from config file and reference transcriptome read from config file
# NOTE: the fastq directory as set in the config is interpreted to be relative to the gExcite working directory
# output: cellranger files
rule cellranger_count_gex:
    input:
        fastqs_dir=WORKDIR
        + "/"
        + config["inputOutput"]["input_fastqs_gex"]
        + "{sample}/",
        reference=config["resources"]["reference_transcriptome"],
    output:
        features_file="results/cellranger_gex/{sample}.features.tsv",
        matrix_file="results/cellranger_gex/{sample}.matrix.mtx",
        barcodes_file="results/cellranger_gex/{sample}.barcodes.tsv",
        web_file=report(
            "results/cellranger_gex/{sample}.web_summary.html",
            caption="../report/cellranger.rst",
            category="preprocessing QC",
        ),
        metrics_file="results/cellranger_gex/{sample}.metrics_summary.csv",
    params:
        cr_out="results/cellranger_gex/",
        variousParams=config["tools"]["cellranger_count_gex"]["variousParams"],
        targetCells=getTargetCellsCellranger_simple,
        mySample="{sample}",
    threads: config["computingResources"]["threads"]["high"]
    log:
        "logs/cellranger_count_gex/{sample}.log",
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["high"],
        runtime=config["computingResources"]["runtime"]["high"],
    benchmark:
        "results/cellranger_gex/benchmark/{sample}.cellranger_count_gex.benchmark"
    # NOTE: cellranger count function cannot specify the output directory, the output is the path you call it from.
    # Therefore, a subshell is used here.
    # Also, unzip and symlink output files in preparation for downstream steps
    shell:
        "(cd {params.cr_out}; "
        "{config[tools][cellranger_count_gex][call]} count "
        "--id={params.mySample} "
        "--sample={params.mySample} "
        "--transcriptome={input.reference} "
        "--fastqs={input.fastqs_dir} "
        "--nosecondary "
        "--localcores={threads} "
        "{params.variousParams}) "
        "&> {log} ; "
        "gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/features.tsv.gz ; "
        "gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ; "
        "gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz ; "
        "ln -frs '{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/features.tsv' '{output.features_file}' ; "
        "ln -frs '{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/matrix.mtx' '{output.matrix_file}' ; "
        "ln -frs '{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/barcodes.tsv' '{output.barcodes_file}' ; "
        "ln -frs '{params.cr_out}{params.mySample}/outs/web_summary.html' '{output.web_file}' ; "
        "ln -frs '{params.cr_out}{params.mySample}/outs/metrics_summary.csv' '{output.metrics_file}' "
