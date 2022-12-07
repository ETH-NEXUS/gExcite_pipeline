
# Import all the basic rules needed from scampi for the analysis one by one.
# Adapt properties where necessary, currently only necessary in create_hdf5.
# All resources and threads are adapted to avoid having a dublicated computingResources section in the config.


# Modify hdf5 rule in order to fit to the directory naming here. 
use rule create_hdf5 from scampi as scampi_create_hdf5 with:
    input:
        genes_file = 'results/cellranger_gex/{sample}.features.tsv',
        matrix_file = 'results/cellranger_gex/{sample}.matrix.mtx',
        barcodes_file = 'results/cellranger_gex/{sample}.barcodes.tsv'
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule identify_doublets from scampi as scampi_identify_doublets with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule filter_genes_and_cells from scampi as scampi_filter_genes_and_cells with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule sctransform_preprocessing from scampi as scampi_sctransform_preprocessing with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule phenograph from scampi as scampi_phenograph with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule prepare_celltyping from scampi as scampi_prepare_celltyping with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule celltyping from scampi as scampi_celltyping with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule remove_atypical_cells from scampi as scampi_remove_atypical with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule gsva from scampi as scampi_gsva with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule plotting from scampi as scampi_plotting with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule gene_exp from scampi as scampi_gene_exp with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule generate_qc_plots_raw from scampi as scampi_generate_qc_plots_raw with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']

use rule generate_cell_type_boxplot from scampi as scampi_generate_cell_type_boxplot with:
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        runtime = config['computingResources']['runtime']['medium']
    threads:
        threads = config['computingResources']['threads']['medium']
