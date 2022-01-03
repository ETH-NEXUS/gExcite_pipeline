# Rules required for the basic cellranger run of sc ADT data
# Lourdes Rosano, January 2020



# Preprocessing to generate library file for ADT cellranger from the provided input parameters
# Library file has a fixed format: fastqs,sample,library_type (with header)
rule generate_library_adt:
    input:
        fastqs_dir=config["inputOutput"]["input_fastqs_adt"] + "{sample}/",
    output:
        library_file="results/pooled_sample/cellranger_adt/{sample}.adt_library.txt",
    params:
        seqRunName=getSeqRunName,
    threads: config["tools"]["generate_library_adt"]["threads"]
    resources: 
        scratch=config["tools"]["generate_library_adt"]["scratch"],
        mem=config["tools"]["generate_library_adt"]["mem"],
        time=config["tools"]["generate_library_adt"]["time"],
        lsfoutfile="results/pooled_sample/cellranger_adt/{sample}.generate_library_adt.lsfout.log",
        lsferrfile="results/pooled_sample/cellranger_adt/{sample}.generate_library_adt.lsferr.log"
    benchmark:
        "results/pooled_sample/cellranger_adt/{sample}.generate_library_adt.benchmark"
    shell:
        'echo -e "fastqs,sample,library_type\n{input.fastqs_dir},{params.seqRunName},Antibody Capture" > {output.library_file}'


# Cellranger call to process the raw ADT samples
rule cellranger_count_adt:
    input:
        library="results/pooled_sample/cellranger_adt/{sample}.adt_library.txt",
        reference=config["resources"]["reference_transcriptome"],
        features_ref=getFeatRefFile,
    output:
        zipped_features_file="results/pooled_sample/cellranger_adt/"
        + "{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        zipped_barcodes_file="results/pooled_sample/cellranger_adt/"
        + "{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        zipped_matrix_file="results/pooled_sample/cellranger_adt/"
        + "{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        features_file="results/pooled_sample/cellranger_adt/{sample}.features.tsv",
        matrix_file="results/pooled_sample/cellranger_adt/{sample}.matrix.mtx",
        barcodes_file="results/pooled_sample/cellranger_adt/{sample}.barcodes.tsv",
        web_file=report(
            "results/pooled_sample/cellranger_adt/{sample}.web_summary.html",
            caption="../report/cellranger_adt.rst",
            category="preprocessing QC",
        ),
        metrics_file="results/pooled_sample/cellranger_adt/{sample}.metrics_summary.csv",
    params:
        cr_out="results/pooled_sample/cellranger_adt/",
        local_cores=config["tools"]["cellranger_count_adt"]["local_cores"],
        variousParams=config["tools"]["cellranger_count_adt"]["variousParams"],
        targetCells=getTargetCellsCellranger,
    threads: config["tools"]["cellranger_count_adt"]["threads"]
    resources:
        lsfoutfile="results/pooled_sample/cellranger_adt/{sample}.cellranger_count_adt.lsfout.log",
        lsferrfile="results/pooled_sample/cellranger_adt/{sample}.cellranger_count_adt.lsferr.log",
        scratch=config["tools"]["cellranger_count_adt"]["scratch"],
        mem=config["tools"]["cellranger_count_adt"]["mem"],
        time=config["tools"]["cellranger_count_adt"]["time"],
    benchmark:
        "results/pooled_sample/cellranger_adt/{sample}.cellranger_count_adt.benchmark"
    # NOTE: cellranger count function cannot specify the output directory, the output it the path you call it from.
    # Therefore, a subshell is used here.
    # Also, unzip and symlink output files in preparation for rule 'create_symlink_adt' or 'create_symlink_adt_nonHashed'
    shell:
        '(cd {params.cr_out}; rm -r {wildcards.sample}/ ; {config[tools][cellranger_count_adt][call]} count --id={wildcards.sample} --transcriptome={input.reference} --localcores={params.local_cores} --libraries={input.library} --feature-ref={input.features_ref} --nosecondary {params.variousParams} {params.targetCells}) && gunzip -c {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/features.tsv.gz > {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/features.tsv ; gunzip -c {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv ; gunzip -c {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz > {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/matrix.mtx && ln -s "{params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/features.tsv" "{output.features_file}"; ln -s "{params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/matrix.mtx" "{output.matrix_file}"; ln -s "{params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv" "{output.barcodes_file}" ; ln -s "{params.cr_out}{wildcards.sample}/outs/web_summary.html" "{output.web_file}" ; ln -s "{params.cr_out}{wildcards.sample}/outs/metrics_summary.csv" "{output.metrics_file}"'

