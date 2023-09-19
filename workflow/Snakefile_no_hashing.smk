# Define minimum Snakemake Version
from snakemake.utils import min_version

min_version("6.12.1")


# Include Config file
configfile: "config/config.yaml"


# Include report functionality
report: "../report/workflow.rst"


# This file includes common functions used in the pipeline
include: "rules/misc_snake.smk"


## Include scampi as a module
module scampi:
    skip_validation:
        True
    snakefile:
        github(
            "ETH-NEXUS/scAmpi_single_cell_RNA",
            path="workflow/snakefile_basic.smk",
            tag="v2.0.7",
        )
    config:
        config["scampi"]


## Include the preprocessing rules
include: "rules/gex_cellranger_no_hashing.smk"
include: "rules/adt_cellranger_no_hashing.smk"
## Include scampi
include: "rules/scampi_module.smk"
## Include citseq rules
include: "rules/adt_analyse_citeseq.smk"


# This rule defines which files should be created
localrules:
    gExcite,


rule gExcite:
    input:
        # Final files from preprocessing
        expand(
            "results/cellranger_gex/{sample}.matrix.mtx",
            sample=getSimpleSampleNames(),
        ),
        expand(
            "results/cellranger_gex/{sample}.barcodes.tsv",
            sample=getSimpleSampleNames(),
        ),
        expand(
            "results/cellranger_gex/{sample}.features.tsv",
            sample=getSimpleSampleNames(),
        ),
        expand(
            "results/cellranger_adt/{sample}.matrix.mtx",
            sample=getSimpleSampleNames(),
        ),
        # List of final files from scampi
        expand(
            "results/counts_raw/{sample}.h5",
            sample=getSimpleSampleNames(),
        ),
        expand(
            "results/counts_filtered/{sample}.doublet_barcodes.txt",
            sample=getSimpleSampleNames(),
        ),
        expand(
            "results/counts_raw/{sample}.raw.histogram_library_sizes.png",
            sample=getSimpleSampleNames(),
        ),
        expand(
            "results/counts_corrected/{sample}.corrected.RDS",
            sample=getSimpleSampleNames(),
        ),
        expand(
            "results/clustering/{sample}.clusters_phenograph.csv",
            sample=getSimpleSampleNames(),
        ),
        expand(
            "results/gene_exp/{sample}.gene_expression_per_cluster.tsv",
            sample=getSimpleSampleNames(),
        ),
        expand(
            "results/plotting/{sample}.celltype_barplot.png",
            sample=getSimpleSampleNames(),
        ),
        expand(
            "results/gsva/{sample}.gsetscore_hm.png",
            sample=getSimpleSampleNames(),
        ),
        # List of final files from citeseq analysis
        expand(
            "results/citeseq_analysis/{sample}/{sample}.GEX_cellrangerADT_SCE.RDS",
            sample=getSimpleSampleNames(),
        ),
    output:
        "results/complete_gExcite.txt",
    resources:
        mem_mb=1000,
        time=1,
    benchmark:
        "results/complete_gExcite.benchmark"
    shell:
        "date > {output} "
