import os, glob, sys, datetime

# Include report functionality
report: "../report/workflow.rst"
# This function adapts the config object to include full path information
include: "./rules/misc_snake.smk"

# input files and folders
SAMPLEMAPPING = config['inputOutput']['sample_map']
INPUTDIR_ADT = config['inputOutput']['input_fastqs_adt']
ROOTDIR = config['inputOutput']['analysis_root_dir']
TMPDIR = config['inputOutput']['analysis_temp_dir']

# Check if the uses specified the proper input and output directories
if not 'INPUTDIR_ADT' in globals():
    print('You have to specify the root directory of the ADT fastq files!')
    sys.exit(1)
if not 'ROOTDIR' in globals():
    print('You have to specify the master root directory where the analysis results will be placed!')
    sys.exit(1)
if not 'TMPDIR' in globals():
    print('You have to specify the root directory where temporary files will be stored!')
    sys.exit(1)

## Create a default outdir for the preprocessing output
OUTDIR = ROOTDIR + 'preprocessing/'

CELLRANGER_ADT_IN = INPUTDIR_ADT
CELLRANGER_ADT_OUT = OUTDIR + 'cellranger_run_adt/'

## Include the rules
include: "./rules/adt_cellranger.smk"
include: "./rules/adt_hashing.smk"

# Require all final files, from both hashed and non-hashed samples (if any)
def getInputFiles(wildcards):
    allFiles = []

    # Cellranger output for hashed samples (in preprocessing folder)
    cr_hashed_adt = expand(CELLRANGER_ADT_OUT + '{sample}.features.tsv', sample = getHashedSampleNames())

    # Cellranger output for non-hashed samples (also in preprocessing folder)
    cr_nonhashed_adt = expand(CELLRANGER_ADT_OUT + '{sample}.features.tsv', sample = getNonHashedSampleNames())

    ####
    # Hashing framework: symlinks to cellranger output (in pooled analysis folders)
    ####
    root_hashed_adt = expand(ROOTDIR + '{sample}/analysis_pooledSample/cellranger_adt/{sample}.features.tsv', sample = getHashedSampleNames())
    # Citeseq output ( in pooled analysis folder)
    citeseq = expand(ROOTDIR + '{sample}/analysis_pooledSample/citeseq_count/{sample}.run_report.yaml',sample = getHashedSampleNames())
    # Hashing output ( in pooled analysis folder)
    hashing = expand(ROOTDIR + '{sample}/analysis_pooledSample/hashing_analysis/{sample}.complete_hashing.txt', sample = getHashedSampleNames())
    ####
    # Non-hashing framework: symlinks to analysis output (in single-sample folder, no pooled)
    ####
    root_nonhashed_adt = expand(ROOTDIR + '{sample}/analysis_{sample}/analysis/cellranger_run_adt/{sample}.features.tsv', sample = getNonHashedSampleNames())

    # Add all strings to list of required files. Empty will be skipped
    # Cellranger results
    for f in cr_hashed_adt:
        allFiles.append(f)
    for f in cr_nonhashed_adt:
        allFiles.append(f)
    # Symlinks to analysis folder
    for f in root_hashed_adt:
        allFiles.append(f)
    for f in root_nonhashed_adt:
        allFiles.append(f)
    # Citeseq
    for f in citeseq:
        allFiles.append(f)
    # Hashing
    for f in hashing:
        allFiles.append(f)

    return allFiles

# This rule defines which files should be created
localrules: scPreprocessing_adt
rule scPreprocessing_adt:
    input:
        inputFiles = getInputFiles
    output:
        OUTDIR + 'complete_sc_preprocessing_adt.txt'
    params:
        lsfoutfile = OUTDIR + 'complete_sc_preprocessing_adt.lsfout.log',
        lsferrfile = OUTDIR + 'complete_sc_preprocessing_adt.lsferr.log',
        mem = '1000',
        scratch = '1000',
        time = '1',
        filter_cellranger = config['tools']['link_filter_cellranger']['call'],
        samplemapFile = SAMPLEMAPPING,
        rootdir = ROOTDIR
    benchmark:
        OUTDIR + 'complete_sc_preprocessing_adt.txt.benchmark'
    shell:
        '{params.filter_cellranger} --sampleMap {params.samplemapFile} --rootdir {params.rootdir} --technology ADT && date > {output}'
