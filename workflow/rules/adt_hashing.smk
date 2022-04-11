# Rules to demultiplex & analyse hashed data. 

from os import listdir
from os.path import isfile, join


# Create Tag File that matches input requirement of CiteSeq Count

rule create_tag_file:
    input:
        tags = getTagFileHashedSamples
    output:
        tagFile = 'results/pooled_samples/citeseq_count/{sample}.tags.tsv'
    params:
        outdir = 'results/pooled_samples/citeseq_count/'
    threads:
        config['computingResources']['threads']['low']
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low']
    log:
        'logs/create_tag_file/{sample}.log'
    benchmark:
        'results/pooled_samples/citeseq_count/benchmark/{sample}.create_tag_file.benchmark'
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
        tags = 'results/pooled_samples/citeseq_count/{sample}.tags.tsv'
    output:
        run_report = report("results/pooled_samples/citeseq_count/{sample}.run_report.yaml", caption="workflow/report/citeseq.rst", category="preprocessing QC")
    params:
        R1 = lambda wildcards: ",".join(list_fastqs(config["inputOutput"]["input_fastqs_adt"]+'{sample}/'.format(sample=wildcards.sample),"R1")),
        R2 = lambda wildcards: ",".join(list_fastqs(config["inputOutput"]["input_fastqs_adt"]+'{sample}/'.format(sample=wildcards.sample),"R2")),
        outdir = 'results/pooled_samples/citeseq_count/{sample}/',
        outfile = 'results/pooled_samples/citeseq_count/{sample}/run_report.yaml',
        variousParams = config['tools']['run_citeseq_count']['variousParams'],
        targetCells = getTargetCellsCiteseqCount
    conda:
        "../envs/run_citeseq_count.yaml"
    resources:
        mem_mb = config['computingResources']['mem']['high'],
        time_min = config['computingResources']['time']['high']
    log:
        'logs/CITE-seq-Count/{sample}.log'
    threads:
        config['computingResources']['threads']['high']
    benchmark:
        'results/pooled_samples/citeseq_count/benchmark/{sample}.run_citeseq_count.benchmark'
    shell:
        'CITE-seq-Count -R1 {params.R1} -R2 {params.R2} {params.variousParams} -o {params.outdir} {params.targetCells} -T {threads} -t {input.tags} ; ln -frs {params.outfile} {output.run_report} &> {log}'

# Hashing Analysis Script requires the zipped cellranger files as input.
# Input: Cellranger output files
# Output: zipped feature files

rule gzip_files_hashingInput:
    input:
        features_file_tmp = 'results/pooled_samples/cellranger_adt/{sample}.features.tsv',
        matrix_file_tmp = 'results/pooled_samples/cellranger_adt/{sample}.matrix.mtx',
        barcodes_file_tmp = 'results/pooled_samples/cellranger_adt/{sample}.barcodes.tsv'
    output:
        features_file = 'results/pooled_samples/cellranger_adt/{sample}_zipped_files/features.tsv.gz',
        matrix_file = 'results/pooled_samples/cellranger_adt/{sample}_zipped_files/matrix.mtx.gz',
        barcodes_file = 'results/pooled_samples/cellranger_adt/{sample}_zipped_files/barcodes.tsv.gz'
    params:
        root_out = 'results/pooled_samples/cellranger_adt/{sample}_zipped_files/'
    resources:
        mem_mb = config['computingResources']['mem']['high'],
        time_min = config['computingResources']['time']['high']
    threads:
        config['computingResources']['threads']['high']
    log:
        'logs/gzip_files_hashingInput/{sample}.log'
    benchmark:
        'results/pooled_samples/cellranger_adt/benchmark/{sample}.create_symlink_adt.benchmark'
    shell:
        'mkdir -p {params.root_out} ; gzip -c "{input.features_file_tmp}" > "{output.features_file}"; gzip -c "{input.matrix_file_tmp}" > "{output.matrix_file}"; gzip -c "{input.barcodes_file_tmp}" > "{output.barcodes_file}"'

# Analyse Hashing: Produces one list / tag containing the barcodes & some qc plots.
# Input: Output of Citeseq and Cellranger ADT run.
# Output: One file / tag + one file / Negatives listing the barcodes found for the tag.

rule Rscript_analyseHashing:
    input:
        citeseq = 'results/pooled_samples/citeseq_count/{sample}.run_report.yaml',
        adt_features = 'results/pooled_samples/cellranger_adt/{sample}_zipped_files/features.tsv.gz',
        adt_matrix = 'results/pooled_samples/cellranger_adt/{sample}_zipped_files/matrix.mtx.gz',
        adt_barcodes = 'results/pooled_samples/cellranger_adt/{sample}_zipped_files/barcodes.tsv.gz',
        tags = getTagFileHashedSamples
    output:
        successFile = 'results/pooled_samples/hashing_analysis/{sample}.complete_hashing.txt'
    conda:
        "../envs/analyse_hashing.yaml"
    params:
        adt_folder = 'results/pooled_samples/cellranger_adt/{sample}_zipped_files/',
        output_prefix = 'results/pooled_samples/hashing_analysis/{sample}',
        quantileThreshold = config['tools']['analyse_hashing']['quantileThreshold'],
        normalisation = config['tools']['analyse_hashing']['normalisation'],
        normalisation_downstream = config['tools']['analyse_hashing']['normalisation_downstream'],
        save_negatives = config['tools']['analyse_hashing']['save_negatives'],
        citeseq_folder = 'results/pooled_samples/citeseq_count/{sample}/umi_count/'
    log:
        'logs/Rscript_analyseHashing/{sample}.log'
    resources:
        mem_mb = config['computingResources']['mem']['high'],
        time_min = config['computingResources']['time']['high']
    threads:
        config['computingResources']['threads']['high']
    benchmark:
        'results/pooled_samples/hashing_analysis/benchmark/{sample}.analyse_hashing.benchmark'
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

rule Rscript_demultiplex_count_matrix:
    input:
        hashingSuccessFile = 'results/pooled_samples/hashing_analysis/{sample}.complete_hashing.txt',
        features_file="results/pooled_samples/cellranger_gex/{sample}.features.tsv",
        matrix_file="results/pooled_samples/cellranger_gex/{sample}.matrix.mtx",
        barcodes_file="results/pooled_samples/cellranger_gex/{sample}.barcodes.tsv"
    output:
        cellranger_gex_mtx = 'results/cellranger_gex/{sample}.{demultiplexed}.matrix.mtx',
        cellranger_gex_feature = 'results/cellranger_gex/{sample}.{demultiplexed}.features.tsv',
        cellranger_gex_barcode = 'results/cellranger_gex/{sample}.{demultiplexed}.barcodes.tsv',
        cellranger_adt_mtx = 'results/cellranger_adt/{sample}.{demultiplexed}.matrix.mtx',
        cellranger_adt_feature = 'results/cellranger_adt/{sample}.{demultiplexed}.features.tsv',
        cellranger_adt_barcode = 'results/cellranger_adt/{sample}.{demultiplexed}.barcodes.tsv',
    conda:
        "../envs/demultiplex_count_matrix.yaml"
    params:
        samplemapFile=config["inputOutput"]["sample_map"]
    log:
        'logs/Rscript_demultiplex_count_matrix/{sample}.{demultiplexed}.log'
    benchmark:
        'results/cellranger_gex/benchmark/{sample}.{demultiplexed}.demultiplex_count_matrix.benchmark'
    resources:
        mem_mb = config['computingResources']['mem']['medium'],
        time_min = config['computingResources']['time']['medium']
    threads:config['computingResources']['threads']['medium']
    shell:
        "Rscript workflow/scripts/demultiplex_count_matrix.R --sampleMap {params.samplemapFile} --sample {wildcards.sample} --rootdir './results/' &> {log}"
