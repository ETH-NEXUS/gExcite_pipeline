
# Import all the basic rules needed from scampi for the analysis one by one.
# Adapt properties where necessary.
# All resources and threads are adapted to avoid having a dublicated computingResources section in the config.


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

use rule celltyping from scampi as scampi_celltyping with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']

use rule remove_atypical_cells from scampi as scampi_remove_atypical with:
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

use rule cell_percent_in_cluster from scampi as scampi_cell_percent_in_cluster with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']


use rule diff_exp_analysis from scampi as scampi_diff_exp_analysis with:
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        threads = config['computingResources']['mediumRequirements']['threads']
