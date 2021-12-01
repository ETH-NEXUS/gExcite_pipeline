# Rules required for the basic cellranger run of sc ADT data
# Lourdes Rosano, January 2020

if not "CELLRANGER_ADT_IN" in globals():
    CELLRANGER_ADT_IN = INPUTDIR_ADT
if not "CELLRANGER_ADT_OUT" in globals():
    CELLRANGER_ADT_OUT = OUTDIR + "cellranger_run_adt/"


# Preprocessing to generate library file for ADT cellranger from the provided input parameters
# Library file has a fixed format: fastqs,sample,library_type (with header)
rule generate_library_adt:
    input:
        fastqs_dir=CELLRANGER_ADT_IN + "{sample}/",
    output:
        library_file=CELLRANGER_ADT_OUT + "{sample}.adt_library.txt",
    params:
        lsfoutfile=CELLRANGER_ADT_OUT + "{sample}.generate_library_adt.lsfout.log",
        lsferrfile=CELLRANGER_ADT_OUT + "{sample}.generate_library_adt.lsferr.log",
        scratch=config["tools"]["generate_library_adt"]["scratch"],
        mem=config["tools"]["generate_library_adt"]["mem"],
        time=config["tools"]["generate_library_adt"]["time"],
        seqRunName=getSeqRunName,
    threads: config["tools"]["generate_library_adt"]["threads"]
    benchmark:
        CELLRANGER_ADT_OUT + "{sample}.generate_library_adt.benchmark"
    shell:
        'echo -e "fastqs,sample,library_type\n{input.fastqs_dir},{params.seqRunName},Antibody Capture" > {output.library_file}'


# Cellranger call to process the raw ADT samples
rule cellranger_count_adt:
    input:
        library=CELLRANGER_ADT_OUT + "{sample}.adt_library.txt",
        reference=config["resources"]["reference_transcriptome"],
        features_ref=getFeatRefFile,
    output:
        zipped_features_file=CELLRANGER_ADT_OUT
        + "{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        zipped_barcodes_file=CELLRANGER_ADT_OUT
        + "{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        zipped_matrix_file=CELLRANGER_ADT_OUT
        + "{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        features_file=CELLRANGER_ADT_OUT + "{sample}.features.tsv",
        matrix_file=CELLRANGER_ADT_OUT + "{sample}.matrix.mtx",
        barcodes_file=CELLRANGER_ADT_OUT + "{sample}.barcodes.tsv",
        web_file=report(
            CELLRANGER_ADT_OUT + "{sample}.web_summary.html",
            caption="../report/cellranger_adt.rst",
            category="preprocessing QC",
        ),
        metrics_file=CELLRANGER_ADT_OUT + "{sample}.metrics_summary.csv",
    params:
        cr_out=CELLRANGER_ADT_OUT,
        local_cores=config["tools"]["cellranger_count_adt"]["local_cores"],
        lsfoutfile=CELLRANGER_ADT_OUT + "{sample}.cellranger_count_adt.lsfout.log",
        lsferrfile=CELLRANGER_ADT_OUT + "{sample}.cellranger_count_adt.lsferr.log",
        scratch=config["tools"]["cellranger_count_adt"]["scratch"],
        mem=config["tools"]["cellranger_count_adt"]["mem"],
        time=config["tools"]["cellranger_count_adt"]["time"],
        variousParams=config["tools"]["cellranger_count_adt"]["variousParams"],
        targetCells=getTargetCellsCellranger,
    threads: config["tools"]["cellranger_count_adt"]["threads"]
    benchmark:
        CELLRANGER_ADT_OUT + "{sample}.cellranger_count_adt.benchmark"
    # NOTE: cellranger count function cannot specify the output directory, the output it the path you call it from.
    # Therefore, a subshell is used here.
    # Also, unzip and symlink output files in preparation for rule 'create_symlink_adt' or 'create_symlink_adt_nonHashed'
    shell:
        '(cd {params.cr_out}; rm -r {wildcards.sample}/ ; {config[tools][cellranger_count_adt][call]} count --id={wildcards.sample} --transcriptome={input.reference} --localcores={params.local_cores} --libraries={input.library} --feature-ref={input.features_ref} --nosecondary {params.variousParams} {params.targetCells}) && gunzip -c {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/features.tsv.gz > {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/features.tsv ; gunzip -c {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv ; gunzip -c {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz > {params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/matrix.mtx && ln -s "{params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/features.tsv" "{output.features_file}"; ln -s "{params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/matrix.mtx" "{output.matrix_file}"; ln -s "{params.cr_out}{wildcards.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv" "{output.barcodes_file}" ; ln -s "{params.cr_out}{wildcards.sample}/outs/web_summary.html" "{output.web_file}" ; ln -s "{params.cr_out}{wildcards.sample}/outs/metrics_summary.csv" "{output.metrics_file}"'


# Hashing framework: create symlinks from preprocessing folder to analysis folder, to be able to continue with downstream analyses
# This rule creates a symlink of a symlink
rule create_symlink_adt:
    input:
        features_file_tmp=CELLRANGER_ADT_OUT + "{sample}.features.tsv",
        matrix_file_tmp=CELLRANGER_ADT_OUT + "{sample}.matrix.mtx",
        barcodes_file_tmp=CELLRANGER_ADT_OUT + "{sample}.barcodes.tsv",
        web_summary_tmp=CELLRANGER_ADT_OUT + "{sample}.web_summary.html",
        metrics_summary_tmp=CELLRANGER_ADT_OUT + "{sample}.metrics_summary.csv",
    output:
        # construct path to corresponding ADT sample in analysis directory
        # NOTE: here, a fixed structure is assumed for the analysis directory
        features_file=ROOTDIR
        + "{sample}/analysis_pooledSample/cellranger_adt/{sample}.features.tsv",
        matrix_file=ROOTDIR
        + "{sample}/analysis_pooledSample/cellranger_adt/{sample}.matrix.mtx",
        barcodes_file=ROOTDIR
        + "{sample}/analysis_pooledSample/cellranger_adt/{sample}.barcodes.tsv",
        web_summary=ROOTDIR
        + "{sample}/analysis_pooledSample/cellranger_adt/{sample}.web_summary.html",
        metrics_summary=ROOTDIR
        + "{sample}/analysis_pooledSample/cellranger_adt/{sample}.metrics_summary.csv",
    params:
        root_out=ROOTDIR + "{sample}/analysis_pooledSample/cellranger_adt/",
        lsfoutfile=ROOTDIR
        + "{sample}/analysis_pooledSample/cellranger_adt/{sample}.create_symlink_adt.lsfout.log",
        lsferrfile=ROOTDIR
        + "{sample}/analysis_pooledSample/cellranger_adt/{sample}.create_symlink_adt.lsferr.log",
        scratch=config["tools"]["create_symlink"]["scratch"],
        mem=config["tools"]["create_symlink"]["mem"],
        time=config["tools"]["create_symlink"]["time"],
    threads: config["tools"]["create_symlink"]["threads"]
    benchmark:
        (
            ROOTDIR
            + "{sample}/analysis_pooledSample/cellranger_adt/{sample}.create_symlink_adt.benchmark"
        )
    shell:
        'mkdir -p {params.root_out} ; ln -s "{input.features_file_tmp}" "{output.features_file}"; ln -s "{input.matrix_file_tmp}" "{output.matrix_file}"; ln -s "{input.barcodes_file_tmp}" "{output.barcodes_file}" ; ln -s "{input.web_summary_tmp}" "{output.web_summary}" ; ln -s "{input.metrics_summary_tmp}" "{output.metrics_summary}"'


# Non-hashing framework: create symlinks from preprocessing folder to final sample analysis folder (keep sample set structure), to be able to continue with downstream analyses
# This rule creates a symlink of a symlink
rule create_symlink_adt_nonHashed:
    input:
        features_file_tmp=CELLRANGER_ADT_OUT + "{sample}.features.tsv",
        matrix_file_tmp=CELLRANGER_ADT_OUT + "{sample}.matrix.mtx",
        barcodes_file_tmp=CELLRANGER_ADT_OUT + "{sample}.barcodes.tsv",
        web_summary_tmp=CELLRANGER_ADT_OUT + "{sample}.web_summary.html",
        metrics_summary_tmp=CELLRANGER_ADT_OUT + "{sample}.metrics_summary.csv",
        zipped_features_file=CELLRANGER_ADT_OUT
        + "{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        zipped_barcodes_file=CELLRANGER_ADT_OUT
        + "{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        zipped_matrix_file=CELLRANGER_ADT_OUT
        + "{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
    output:
        # construct path to corresponding sample-specific ADT analysis directory
        # NOTE: here, a fixed structure is assumed for the analysis directory
        features_file=ROOTDIR
        + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}.features.tsv",
        matrix_file=ROOTDIR
        + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}.matrix.mtx",
        barcodes_file=ROOTDIR
        + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}.barcodes.tsv",
        web_summary=ROOTDIR
        + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}.web_summary.html",
        metrics_summary=ROOTDIR
        + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}.metrics_summary.csv",
        zipped_features_file=ROOTDIR
        + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        zipped_matrix_file=ROOTDIR
        + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}/outs/filtered_feature_bc_matrix/matrix.tsv.gz",
        zipped_barcodes_file=ROOTDIR
        + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    params:
        root_out=ROOTDIR + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/",
        lsfoutfile=ROOTDIR
        + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}.create_symlink_adt_nonHashed.lsfout.log",
        lsferrfile=ROOTDIR
        + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}.create_symlink_adt_nonHashed.lsferr.log",
        scratch=config["tools"]["create_symlink"]["scratch"],
        mem=config["tools"]["create_symlink"]["mem"],
        time=config["tools"]["create_symlink"]["time"],
    threads: config["tools"]["create_symlink"]["threads"]
    benchmark:
        (
            ROOTDIR
            + "{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}.create_symlink_adt_nonHashed.benchmark"
        )
    shell:
        'mkdir -p {params.root_out} ; ln -s "{input.features_file_tmp}" "{output.features_file}"; ln -s "{input.matrix_file_tmp}" "{output.matrix_file}"; ln -s "{input.barcodes_file_tmp}" "{output.barcodes_file}" ; ln -s "{input.web_summary_tmp}" "{output.web_summary}" ; ln -s "{input.metrics_summary_tmp}" "{output.metrics_summary}" ; ln -s "{input.zipped_features_file}" "{output.zipped_features_file}"; ln -s "{input.zipped_barcodes_file}" "{output.zipped_barcodes_file}" ;ln -s "{input.zipped_matrix_file}" "{output.zipped_matrix_file}"'
