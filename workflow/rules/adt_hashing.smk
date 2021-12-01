# Rules to demultiplex & analyse hashed data. 

from os import listdir
from os.path import isfile, join


def list_fastqs(base,Rvalue):
    """Returns comma separated list of fastqs."""
    onlyfiles = [base + f for f in listdir(base) if isfile(join(base, f)) and Rvalue in f]
    return onlyfiles


# Create Tag File that matches input requirement of CiteSeq Count

rule createTagFile:
    input:
        tags = getTagFileHashedSamples
    output:
        tagFile = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.tags.tsv',
    params:
        lsfoutfile = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.create_tag_file.lsfout.log',
        lsferrfile = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.create_tag_file.lsferr.log',
        scratch = config['tools']['create_tag_file']['scratch'],
        mem = config['tools']['create_tag_file']['mem'],
        time = config['tools']['create_tag_file']['time'],
        outdir = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/'
    threads:
        config['tools']['create_tag_file']['threads']
    benchmark:
        ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.create_tag_file.benchmark'
    run:
        outTag = open(output.tagFile,"w")
        infile = open(input.tags,"r")
        for line in infile:
            outTag.write(line.split(",")[0]+","+line.split(",")[1] + "\n")


# Cite-Seq Count: Counts Hashtag per barcode.
# Input: Comma separated list of fastqs.
# Output: Matrix, barcode & features file

if not 'CITESEQ_IN' in globals():
    CITESEQ_IN = INPUTDIR_ADT

rule run_citeseq_count:
    input:
        R1 = lambda wildcards: list_fastqs(CITESEQ_IN+'{sample}/'.format(sample=wildcards.sample),"R1"),
        R2 = lambda wildcards: list_fastqs(CITESEQ_IN+'{sample}/'.format(sample=wildcards.sample),"R2"),
        tags = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.tags.tsv'
    output:
        run_report = report(ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.run_report.yaml', caption="../report/citeseq.rst", category="preprocessing QC"),
    params:
        lsfoutfile = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.run_citeseq_count.lsfout.log',
        lsferrfile = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.run_citeseq_count.lsferr.log',
        scratch = config['tools']['run_citeseq_count']['scratch'],
        mem = config['tools']['run_citeseq_count']['mem'],
        time = config['tools']['run_citeseq_count']['time'],
        outdir = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/',
        outfile = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/run_report.yaml',
        variousParams = config['tools']['run_citeseq_count']['variousParams'],
        targetCells = getTargetCellsCiteseqCount,
        call = config['tools']['run_citeseq_count']['call']
    threads:
        config['tools']['run_citeseq_count']['threads']
    benchmark:
        ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.run_citeseq_count.benchmark'
    run:
        R1=",".join(input.R1)
        R2=",".join(input.R2)
        shell("{params.call} -R1 {R1} -R2 {R2} {params.variousParams} -o {params.outdir} {params.targetCells} -T {threads} -t {input.tags}; ln -fs {params.outfile} {output.run_report}")

if not 'CELLRANGER_ADT_OUT' in globals():
    CELLRANGER_ADT_OUT = OUTDIR + 'cellranger_run_adt/'

# Hashing Analysis Script requires the zipped cellranger files as input.
# Input: Cellranger output files
# Output: zipped feature files

rule create_symlink_hashing_input:
    input:
        features_file_tmp = CELLRANGER_ADT_OUT + '{sample}.features.tsv',
        matrix_file_tmp = CELLRANGER_ADT_OUT +'{sample}.matrix.mtx',
        barcodes_file_tmp = CELLRANGER_ADT_OUT + '{sample}.barcodes.tsv'
    output:
        # construct path to corresponding ADT sample in analysis directory
        # NOTE: here, a fixed structure is assumed for the analysis directory
        features_file = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/features.tsv.gz',
        matrix_file = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/matrix.mtx.gz',
        barcodes_file = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/barcodes.tsv.gz',
    params:
        root_out = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/',
        lsfoutfile = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}.create_symlink_adt.lsfout.log',
        lsferrfile = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}.create_symlink_adt.lsferr.log',
        scratch = config['tools']['create_symlink']['scratch'],
        mem = config['tools']['create_symlink']['mem'],
        time = config['tools']['create_symlink']['time']
    threads:
        config['tools']['create_symlink']['threads']
    benchmark:
        ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}.create_symlink_adt.benchmark'
    shell:
        'mkdir -p {params.root_out} ; gzip -c "{input.features_file_tmp}" > "{output.features_file}"; gzip -c "{input.matrix_file_tmp}" > "{output.matrix_file}"; gzip -c "{input.barcodes_file_tmp}" > "{output.barcodes_file}"'


if not 'CELLRANGER_NOVA_ADT_OUT' in globals():
    CELLRANGER_NOVA_ADT_OUT = OUTDIR + 'cellranger_run_index_hopping_adt/'

# Hashing Analysis Script requires the zipped cellranger files of novaseq data as input.
# Input: Cellranger output files
# Output: zipped feature files

rule create_symlink_hashing_input_nova:
    input:
        features_file_tmp = CELLRANGER_NOVA_ADT_OUT + '{sample}.features.tsv',
        matrix_file_tmp = CELLRANGER_NOVA_ADT_OUT +'{sample}.matrix.mtx',
        barcodes_file_tmp = CELLRANGER_NOVA_ADT_OUT + '{sample}.barcodes.tsv'
    output:
        # construct path to corresponding ADT sample in analysis directory
        # NOTE: here, a fixed structure is assumed for the analysis directory
        features_file = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/features.tsv.gz',
        matrix_file = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/matrix.mtx.gz',
        barcodes_file = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/barcodes.tsv.gz',
    params:
        root_out = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/',
        lsfoutfile = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}.create_symlink_adt.lsfout.log',
        lsferrfile = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}.create_symlink_adt.lsferr.log',
        scratch = config['tools']['create_symlink']['scratch'],
        mem = config['tools']['create_symlink']['mem'],
        time = config['tools']['create_symlink']['time']
    threads:
        config['tools']['create_symlink']['threads']
    benchmark:
        ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}.create_symlink_adt.benchmark'
    shell:
        'mkdir -p {params.root_out} ; gzip -c "{input.features_file_tmp}" > "{output.features_file}"; gzip -c "{input.matrix_file_tmp}" > "{output.matrix_file}"; gzip -c "{input.barcodes_file_tmp}" > "{output.barcodes_file}"'

# Analyse Hashing: Produces one list / tag containing the barcodes & some qc plots.
# Input: Output of Citeseq and Cellranger ADT run.
# Output: One file / tag + one file / Negatives listing the barcodes found for the tag.

rule analyse_hashing:
    input:
        citeseq = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.run_report.yaml',
        adt_features = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/features.tsv.gz',
        adt_matrix = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/matrix.mtx.gz',
        adt_barcodes = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/barcodes.tsv.gz',
        tags = getTagFileHashedSamples
    output:
        successFile = ROOTDIR + '{sample}/analysis_pooledSample/hashing_analysis/{sample}.complete_hashing.txt'
    params:
        lsfoutfile = ROOTDIR + '{sample}/analysis_pooledSample/hashing_analysis/{sample}.analyse_hashing.lsfout.log',
        lsferrfile = ROOTDIR + '{sample}/analysis_pooledSample/hashing_analysis/{sample}.analyse_hashing.lsferr.log',
        scratch = config['tools']['analyse_hashing']['scratch'],
        mem = config['tools']['analyse_hashing']['mem'],
        time = config['tools']['analyse_hashing']['time'],
        adt_folder = ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}_zipped_files/',
        output_prefix = ROOTDIR + '{sample}/analysis_pooledSample/hashing_analysis/{sample}',
        quantileThreshold = config['tools']['analyse_hashing']['quantileThreshold'],
        normalisation = config['tools']['analyse_hashing']['normalisation'],
        normalisation_downstream = config['tools']['analyse_hashing']['normalisation_downstream'],
        save_negatives = config['tools']['analyse_hashing']['save_negatives'],
        call = config['tools']['analyse_hashing']['call'],
        citeseq_folder = ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/umi_count/'
    threads:
        config['tools']['analyse_hashing']['threads']
    benchmark:
        ROOTDIR + '{sample}/analysis_pooledSample/hashing_analysis/{sample}.analyse_hashing.benchmark'
    shell:
        '{params.call} ' +
        '--citeseq_in {params.citeseq_folder} ' +
        '--adt_barcodes_in {params.adt_folder} ' +
        '--quantile_threshold {params.quantileThreshold} ' +
        '--sampleTagMap {input.tags} ' +
        '--normalisation {params.normalisation} ' +
        '--normalisation_downstream {params.normalisation_downstream} ' +
        '--save_negatives {params.save_negatives} ' +
        '--output_prefix {params.output_prefix} '+
        '&& date > {output.successFile}'
