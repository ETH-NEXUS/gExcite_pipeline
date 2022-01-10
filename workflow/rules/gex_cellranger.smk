# Cellranger call to process the raw GEX samples
rule cellranger_count_gex:
    input:
        fastqs_dir=config["inputOutput"]["input_fastqs_gex"] + "{sample}/",
        reference=config["resources"]["reference_transcriptome"]
    output:
        features_file="results/pooled_sample/cellranger_gex/{sample}.features.tsv",
        matrix_file="results/pooled_sample/cellranger_gex/{sample}.matrix.mtx",
        barcodes_file="results/pooled_sample/cellranger_gex/{sample}.barcodes.tsv",
        web_file=report(
            "results/pooled_sample/cellranger_gex/{sample}.web_summary.html",
            caption="../report/cellranger.rst",
            category="preprocessing QC",
        ),
        metrics_file="results/pooled_sample/cellranger_gex/{sample}.metrics_summary.csv",
    params:
        cr_out="results/pooled_sample/cellranger_gex/",
        local_cores=config["tools"]["cellranger_count_gex"]["local_cores"],
        variousParams=config["tools"]["cellranger_count_gex"]["variousParams"],
        targetCells=getTargetCellsCellranger,
        metrics_summary = "results/pooled_sample/{sample}.metrics_summary.csv",
        web_summary = "results/pooled_sample/{sample}.web_summary.html",
	mySample = '{sample}' 
    threads: config["tools"]["cellranger_count_gex"]["threads"]
    resources:
        mem_mb=config["tools"]["cellranger_count_gex"]["mem"],
        time_min=config["tools"]["cellranger_count_gex"]["time"]
    benchmark:
        "results/pooled_sample/cellranger_gex/{sample}.cellranger_count_gex.benchmark"
    # NOTE: cellranger count function cannot specify the output directory, the output it the path you call it from.
    # Therefore, a subshell is used here.
    # Also, unzip and symlink output files in preparation for downstream steps
    shell:
         '(cd {params.cr_out}; {config[tools][cellranger_count_gex][call]} count --id={params.mySample} --sample={params.mySample} --transcriptome={input.reference} --localcores={params.local_cores} --fastqs={input.fastqs_dir} --nosecondary {params.variousParams}); gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/features.tsv.gz ; gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ; gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz ; ln -s "{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/features.tsv" "{output.features_file}"; ln -s "{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/matrix.mtx" "{output.matrix_file}"; ln -s "{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/barcodes.tsv" "{output.barcodes_file}" ; ln -s "{params.cr_out}{params.mySample}/outs/web_summary.html" "{params.web_summary}" ; ln -s "{params.cr_out}{params.mySample}/outs/metrics_summary.csv" "{params.metrics_summary}"'
