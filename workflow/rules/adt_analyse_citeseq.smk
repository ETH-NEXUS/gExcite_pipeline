#Rule to Analyse ADT Data in combination with GEXdata
rule analyse_citeseq:
    input:
        RDS = 'results/atypical_removed/{sample}.genes_cells_filtered.corrected.atypical_removed.RDS' ,
        CellrangerADT = 'results/cellranger_adt/{sample}.matrix.mtx',
        h5 = 'results/counts_corrected/{sample}.genes_cells_filtered.corrected.variable_genes.h5',
    output:
        RDS = 'results/citeseq_analysis/{sample}.GEX_cellrangerADT_SCE.RDS',
        completeFile = 'results/citeseq_analysis/{sample}.citeseq_analysis.complete.txt' 
    conda:
        "../envs/adt_analyse_citeseq.yaml"
    params:
        ADTFolder = 'results/cellranger_adt/{sample}/outs/',
        colorConfig = config['scampi']['resources']['colour_config'] ,
        lookup = config['resources']['adt_lookup'],
        threshold = config['resources']['adt_thresholds'],
        numberVariableGenes = config['tools']['analyse_citeseq']['numberVariableGenes'],
        outdir = 'results/citeseq_analysis/' 
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:config['computingResources']['mediumRequirements']['threads']
    shell:
        "Rscript ../scripts/analyse_citeseq.R " +
        "--RDS {input.RDS}" +
        "--cellrangerADT {params.ADTFolder}" +
        "--h5 {input.h5}" +
        "--colorConfig {params.colorConfig}" +
        "--lookup {params.lookup}" +
        "--threshold {params.threshold}" +
        "--threads {threads}" +
        "--sampleName {wildcards.sample}" +
        "--number_variable_genes {params.numberVariableGenes} " +
        "--output_directory {params.outdir} && touch {output.completeFile}"
