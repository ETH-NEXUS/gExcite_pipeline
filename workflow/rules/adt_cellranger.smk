# Rules required for the basic cellranger run of sc ADT data

import os

WORKDIR = os.getcwd()


# Preprocessing to generate library file for ADT cellranger from the provided input parameters
# Library file has a fixed format: fastqs,sample,library_type (with header)
rule create_library_file_adt:
    input:
        fastqs_dir=config["inputOutput"]["input_fastqs_adt"] + "{sample_set}/",
    output:
        library_file=WORKDIR
        + "/results/pooled_samples/cellranger_adt/{sample_set}.adt_library.txt",
    params:
        seqRunName=getSeqRunName,
    threads: config["computingResources"]["threads"]["low"]
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["low"],
        runtime=config["computingResources"]["runtime"]["low"],
    log:
        "logs/create_library_file_adt/{sample_set}.log",
    benchmark:
        "results/pooled_samples/cellranger_adt/benchmark/{sample_set}.generate_library_adt.benchmark"
    shell:
        r"""
        echo -e "fastqs,sample,library_type\n{input.fastqs_dir},{wildcards.sample_set},Antibody Capture" > {output.library_file}
        """


# Cellranger call to process the raw ADT samples
rule cellranger_count_adt:
    input:
        library=WORKDIR
        + "/results/pooled_samples/cellranger_adt/{sample_set}.adt_library.txt",
        reference=config["resources"]["reference_transcriptome"],
        features_ref=getFeatRefFile,
    output:
        features_file="results/pooled_samples/cellranger_adt/{sample_set}.features.tsv",
        matrix_file="results/pooled_samples/cellranger_adt/{sample_set}.matrix.mtx",
        barcodes_file="results/pooled_samples/cellranger_adt/{sample_set}.barcodes.tsv",
        web_file=report(
            "results/pooled_samples/cellranger_adt/{sample_set}.web_summary.html",
            caption="../report/cellranger_adt.rst",
            category="preprocessing QC",
        ),
        metrics_file="results/pooled_samples/cellranger_adt/{sample_set}.metrics_summary.csv",
    params:
        cr_out="results/pooled_samples/cellranger_adt/",
        variousParams=config["tools"]["cellranger_count_adt"]["variousParams"],
        targetCells=getTargetCellsCellranger,
        sample="{sample_set}",
    threads: config["computingResources"]["threads"]["high"]
    log:
        "logs/cellranger_count_adt/{sample_set}.log",
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["high"],
        runtime=config["computingResources"]["runtime"]["high"],
    benchmark:
        "results/pooled_samples/cellranger_adt/benchmark/{sample_set}.cellranger_count_adt.benchmark"
    # NOTE: cellranger count function cannot specify the output directory, the output it the path you call it from.
    # Therefore, a subshell is used here.
    shell:
        "(cd {params.cr_out} ; "
        "{config[tools][cellranger_count_adt][call]} count "
        "--id={params.sample} "
        "--transcriptome={input.reference} "
        "--libraries={input.library} "
        "--feature-ref={input.features_ref} "
        "--nosecondary {params.variousParams} {params.targetCells}) "
        "&> {log} "
        "&& gunzip -c {params.cr_out}{params.sample}/outs/filtered_feature_bc_matrix/features.tsv.gz > {params.cr_out}{params.sample}/outs/filtered_feature_bc_matrix/features.tsv ; "
        "gunzip -c {params.cr_out}{params.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > {params.cr_out}{params.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv ; "
        "gunzip -c {params.cr_out}{params.sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz > {params.cr_out}{params.sample}/outs/filtered_feature_bc_matrix/matrix.mtx "
        "&& ln -rs '{params.cr_out}{params.sample}/outs/filtered_feature_bc_matrix/features.tsv' '{output.features_file}' ; "
        "ln -rs '{params.cr_out}{params.sample}/outs/filtered_feature_bc_matrix/matrix.mtx' '{output.matrix_file}' ; "
        "ln -rs '{params.cr_out}{params.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv' '{output.barcodes_file}' ; "
        "ln -rs '{params.cr_out}{params.sample}/outs/web_summary.html' '{output.web_file}' ; "
        "ln -rs '{params.cr_out}{params.sample}/outs/metrics_summary.csv' '{output.metrics_file}' "
