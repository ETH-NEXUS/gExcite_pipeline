
# Import all the basic rules needed from scampi for the analysis one by one.
# Adapt properties where necessary.
# All resources and threads are adapted to avoid having a dublicated computingResources section in the config.


# Cellranger seems to be imported as well adapt it that dry run works. 
use rule cellranger_count from scampi as scampi_cellranger_count with:
    input:
        fastqs_dir = config['inputOutput']['input_fastqs_gex'],
        reference = config['resources']['reference_transcriptome']
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']


# Modify hdf5 rule in order to fit to the directory naming here. 
use rule create_hdf5 from scampi as scampi_create_hdf5 with:
    input:
        genes_file = 'results/cellranger_gex/{sample}.features.tsv',
        matrix_file = 'results/cellranger_gex/{sample}.matrix.mtx',
        barcodes_file = 'results/cellranger_gex/{sample}.barcodes.tsv'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule identify_doublets from scampi as scampi_identify_doublets with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule filter_genes_and_cells from scampi as scampi_filter_genes_and_cells with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule sctransform_preprocessing from scampi as scampi_sctransform_preprocessing with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule phenograph from scampi as scampi_phenograph with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule prepare_celltyping from scampi as scampi_prepare_celltyping with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule cell_type_classification from scampi as scampi_cell_type_classification with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule remove_atypical from scampi as scampi_remove_atypical with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule gsva from scampi as scampi_gsva with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule plotting from scampi as scampi_plotting with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

# Adapted the params here to avoid having a dublicated inputOutput section in the config file
use rule diff_exp_genes from scampi as scampi_diff_exp_genes with:
    params: 
        malignant = config['inputOutput']['malignant_cell_type'],
        sampleName = '{sample}',
        threshold_comparison = config['scampi']['tools']['diff_exp']['threshold_comparison'],
        fdr_cut = config['scampi']['tools']['diff_exp']['fdr_cut'],
        fc_cut = config['scampi']['tools']['diff_exp']['fc_cut'],
        mindiff2second = config['scampi']['tools']['diff_exp']['mindiff2second'],
        minNumberNonMalignant = config['scampi']['tools']['diff_exp']['minNumberNonMalignant'],
        outpath = 'results/diff_exp/'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule gene_exp from scampi as scampi_gene_exp with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule generate_qc_plots from scampi as scampi_generate_qc_plots with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule generate_cell_type_boxplot from scampi as scampi_generate_cell_type_boxplot with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']


