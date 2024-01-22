
import os
WORKDIR = os.getcwd()

# Cellranger call to process the raw GEX samples
# input: fastq directory from config file and reference transcriptome read from config file
# NOTE: the fastq directory as set in the config is interpreted to be relative to the gExcite working directory
# output: cellranger files
rule cellranger_count_gex:
    input:
        fastqs_dir=WORKDIR + "/" + config["inputOutput"]["input_fastqs_gex"] + "{sample_set}/",
        reference=config["resources"]["reference_transcriptome"],
    output:
        features_file="results/pooled_samples/cellranger_gex/{sample_set}.features.tsv",
        matrix_file="results/pooled_samples/cellranger_gex/{sample_set}.matrix.mtx",
        barcodes_file="results/pooled_samples/cellranger_gex/{sample_set}.barcodes.tsv",
        web_file=report(
            "results/pooled_samples/cellranger_gex/{sample_set}.web_summary.html",
            caption="../report/cellranger.rst",
            category="preprocessing QC",
        ),
        metrics_file="results/pooled_samples/cellranger_gex/{sample_set}.metrics_summary.csv",
    params:
        cr_out="results/pooled_samples/cellranger_gex/",
        variousParams=config["tools"]["cellranger_count_gex"]["variousParams"],
        targetCells=getTargetCellsCellranger,
        mySample="{sample_set}",
    threads: config["computingResources"]["threads"]["high"]
    log:
        "logs/cellranger_count_gex/{sample_set}.log",
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["high"],
        runtime=config["computingResources"]["runtime"]["high"],
    benchmark:
        "results/pooled_samples/cellranger_gex/benchmark/{sample_set}.cellranger_count_gex.benchmark"
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
        "gzip -dk {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/features.tsv.gz ; "
        "gzip -dk {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ; "
        "gzip -dk {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz ; "
        "ln -frs '{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/features.tsv' '{output.features_file}' ; "
        "ln -frs '{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/matrix.mtx' '{output.matrix_file}' ; "
        "ln -frs '{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/barcodes.tsv' '{output.barcodes_file}' ; "
        "ln -frs '{params.cr_out}{params.mySample}/outs/web_summary.html' '{output.web_file}' ; "
        "ln -frs '{params.cr_out}{params.mySample}/outs/metrics_summary.csv' '{output.metrics_file}' "
