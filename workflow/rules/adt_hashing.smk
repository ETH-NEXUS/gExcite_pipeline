# Rules to demultiplex & analyse hashed data. 
from os import listdir
from os.path import isfile, join


# Create Tag File that matches input requirement of CiteSeq Count

rule create_tag_file:
    input:
        tags = getTagFileHashedSamples
    output:
        tagFile = 'results/pooled_samples/citeseq_count/{sample_set}.tags.tsv'
    params:
        outdir = 'results/pooled_samples/citeseq_count/'
    threads:
        config['computingResources']['threads']['low']
    resources:
        mem_mb = config['computingResources']['mem_mb']['low'],
        runtime = config['computingResources']['runtime']['low']
    log:
        'logs/create_tag_file/{sample_set}.log'
    benchmark:
        'results/pooled_samples/citeseq_count/benchmark/{sample_set}.create_tag_file.benchmark'
    run:
        outTag = open(output.tagFile,"w")
        infile = open(input.tags,"r")
        for line in infile:
            outTag.write(line.split(",")[0]+","+line.split(",")[1] + "\n")


# Cite-Seq Count: Counts Hashtag per barcode.
# Input: Comma separated list of fastqs.
# Output: Matrix, barcode & features file

rule CITE_seq_Count:
    input:
        tags = 'results/pooled_samples/citeseq_count/{sample_set}.tags.tsv'
    output:
        run_report = report("results/pooled_samples/citeseq_count/{sample_set}.run_report.yaml", caption="workflow/report/citeseq.rst", category="preprocessing QC")
    params:
        R1 = lambda wildcards: ",".join(list_fastqs(config["inputOutput"]["input_fastqs_adt"]+'{sample_set}/'.format(sample_set=wildcards.sample_set),"R1")),
        R2 = lambda wildcards: ",".join(list_fastqs(config["inputOutput"]["input_fastqs_adt"]+'{sample_set}/'.format(sample_set=wildcards.sample_set),"R2")),
        outdir = 'results/pooled_samples/citeseq_count/{sample_set}/',
        outfile = 'results/pooled_samples/citeseq_count/{sample_set}/run_report.yaml',
        variousParams = config['tools']['run_citeseq_count']['variousParams'],
        targetCells = getTargetCellsCiteseqCount
    conda:
        "../envs/run_citeseq_count.yaml"
    resources:
        mem_mb = config['computingResources']['mem_mb']['high'],
        runtime = config['computingResources']['runtime']['high']
    log:
        'logs/CITE-seq-Count/{sample_set}.log'
    threads:
        config['computingResources']['threads']['high']
    benchmark:
        'results/pooled_samples/citeseq_count/benchmark/{sample_set}.run_citeseq_count.benchmark'
    shell:
        'CITE-seq-Count -R1 {params.R1} -R2 {params.R2} {params.variousParams} -o {params.outdir} {params.targetCells} -T {threads} -t {input.tags} ; ln -frs {params.outfile} {output.run_report} &> {log}'

# Hashing Analysis Script requires the zipped cellranger files as input.
# Input: Cellranger output files
# Output: zipped feature files

rule gzip_files_hashingInput:
    input:
        features_file_tmp = 'results/pooled_samples/cellranger_adt/{sample_set}.features.tsv',
        matrix_file_tmp = 'results/pooled_samples/cellranger_adt/{sample_set}.matrix.mtx',
        barcodes_file_tmp = 'results/pooled_samples/cellranger_adt/{sample_set}.barcodes.tsv'
    output:
        features_file = 'results/pooled_samples/cellranger_adt/{sample_set}_zipped_files/features.tsv.gz',
        matrix_file = 'results/pooled_samples/cellranger_adt/{sample_set}_zipped_files/matrix.mtx.gz',
        barcodes_file = 'results/pooled_samples/cellranger_adt/{sample_set}_zipped_files/barcodes.tsv.gz'
    params:
        root_out = 'results/pooled_samples/cellranger_adt/{sample_set}_zipped_files/'
    resources:
        mem_mb = config['computingResources']['mem_mb']['high'],
        runtime = config['computingResources']['runtime']['high']
    threads:
        config['computingResources']['threads']['high']
    log:
        'logs/gzip_files_hashingInput/{sample_set}.log'
    benchmark:
        'results/pooled_samples/cellranger_adt/benchmark/{sample_set}.create_symlink_adt.benchmark'
    shell:
        'mkdir -p {params.root_out} ; gzip -c "{input.features_file_tmp}" > "{output.features_file}"; gzip -c "{input.matrix_file_tmp}" > "{output.matrix_file}"; gzip -c "{input.barcodes_file_tmp}" > "{output.barcodes_file}"'

# Analyse Hashing: Produces one list / tag containing the barcodes & some qc plots.
# Input: Output of Citeseq and Cellranger ADT run.
# Output: One file / tag + one file / Negatives listing the barcodes found for the tag.

rule Rscript_analyseHashing:
    input:
        citeseq = 'results/pooled_samples/citeseq_count/{sample_set}.run_report.yaml',
        adt_features = 'results/pooled_samples/cellranger_adt/{sample_set}_zipped_files/features.tsv.gz',
        adt_matrix = 'results/pooled_samples/cellranger_adt/{sample_set}_zipped_files/matrix.mtx.gz',
        adt_barcodes = 'results/pooled_samples/cellranger_adt/{sample_set}_zipped_files/barcodes.tsv.gz',
        tags = getTagFileHashedSamples
    output:
        successFile = 'results/pooled_samples/hashing_analysis/{sample_set}.complete_hashing.txt',
    conda:
        "../envs/analyse_hashing.yaml"
    params:
        adt_folder = 'results/pooled_samples/cellranger_adt/{sample_set}_zipped_files/',
        output_prefix = 'results/pooled_samples/hashing_analysis/{sample_set}',
        quantileThreshold = config['tools']['analyse_hashing']['quantileThreshold'],
        normalisation = config['tools']['analyse_hashing']['normalisation'],
        normalisation_downstream = config['tools']['analyse_hashing']['normalisation_downstream'],
        save_negatives = config['tools']['analyse_hashing']['save_negatives'],
        citeseq_folder = 'results/pooled_samples/citeseq_count/{sample_set}/umi_count/'
    log:
        'logs/Rscript_analyseHashing/{sample_set}.log'
    resources:
        mem_mb = config['computingResources']['mem_mb']['high'],
        runtime = config['computingResources']['runtime']['high']
    threads:
        config['computingResources']['threads']['high']
    benchmark:
        'results/pooled_samples/hashing_analysis/benchmark/{sample_set}.analyse_hashing.benchmark'
    shell:
        'Rscript ./workflow/scripts/analyseHashing.R ' +
        '--citeseq_in {params.citeseq_folder} ' +
        '--adt_barcodes_in {params.adt_folder} ' +
        '--quantile_threshold {params.quantileThreshold} ' +
        '--sampleTagMap {input.tags} ' +
        '--normalisation {params.normalisation} ' +
        '--normalisation_downstream {params.normalisation_downstream} ' +
        '--save_negatives {params.save_negatives} ' +
        '--output_prefix {params.output_prefix} '+
        '&> {log}; touch {output.successFile}'

# Rule to demultiplex the count matrix
# Input: Hashing success File from the rule Rscript_analyseHashing and the cellranger files
# Output: cellranger gex and adt files for every sample in sampleset (note! Rule runs once / sample)
rule Rscript_demultiplex_count_matrix:
    input:
        successFile = 'results/pooled_samples/hashing_analysis/{sample_set}.complete_hashing.txt',
        features_file="results/pooled_samples/cellranger_gex/{sample_set}.features.tsv",
        matrix_file="results/pooled_samples/cellranger_gex/{sample_set}.matrix.mtx",
        barcodes_file="results/pooled_samples/cellranger_gex/{sample_set}.barcodes.tsv"
    output:
         features_file = 'results/cellranger_gex/{sample_set}.{sample}.features.tsv',
         matrix_file = 'results/cellranger_gex/{sample_set}.{sample}.matrix.mtx',
         barcodes_file = 'results/cellranger_gex/{sample_set}.{sample}.barcodes.tsv',
         features_file_adt = 'results/cellranger_adt/{sample_set}.{sample}.features.tsv',
         matrix_file_adt = 'results/cellranger_adt/{sample_set}.{sample}.matrix.mtx',
         barcodes_file_adt = 'results/cellranger_adt/{sample_set}.{sample}.barcodes.tsv',
    conda:
        "../envs/demultiplex_count_matrix.yaml"
    params:
        samplemapFile=config["inputOutput"]["sample_map"]
    log:
         'logs/Rscript_demultiplex_count_matrix/{sample_set}.{sample}.log'
    benchmark:
         'results/cellranger_gex/benchmark/{sample_set}.{sample}.demultiplex_count_matrix.benchmark'
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:config['computingResources']['threads']['medium']
    shell:
        "Rscript workflow/scripts/demultiplex_count_matrix.R --sampleMap {params.samplemapFile} --sample {wildcards.sample} --sampleSet {wildcards.sample_set} &> {log} "
