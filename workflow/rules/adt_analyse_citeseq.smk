#Rule to create an initial threshold file
rule create_initial_threshold_file:
    input:
        CellrangerADT = 'results/cellranger_adt/{sample}.features.tsv',
    output:
        thresholds = 'results/citeseq_analysis/{sample}.thresholds.tsv',
    resources:
        mem_mb = config['computingResources']['lowRequirements']['mem'],
        time_min = config['computingResources']['lowRequirements']['time']
    threads:config['computingResources']['lowRequirements']['threads']
    shell:
        "echo -e 'Antibody\t{wildcards.sample}' > {output.thresholds} ; " +
	"awk ' {{ print $2 }}' {input.CellrangerADT} | sed 's/$/\t0/' >> {output.thresholds}"
        

#Rule to Analyse ADT Data in combination with GEXdata
rule analyse_citeseq:
    input:
        RDS = 'results/atypical_removed/{sample}.atypical_removed.RDS' ,
        CellrangerADT = 'results/cellranger_adt/{sample}.matrix.mtx',
        thresholds = 'results/citeseq_analysis/{sample}.thresholds.tsv',
        h5 = 'results/counts_corrected/{sample}.corrected.variable_genes.h5',
    output:
        RDS = 'results/citeseq_analysis/{sample}/{sample}.GEX_cellrangerADT_SCE.RDS',
        completeFile = 'results/citeseq_analysis/{sample}/{sample}.citeseq_analysis.complete.txt' 
    conda:
        "../envs/adt_analyse_citeseq.yaml"
    params:
        ADTFolder = 'results/cellranger_adt/{sample}/outs/',
        colorConfig = config['scampi']['resources']['colour_config'] ,
        lookup = config['resources']['adt_lookup'],
        numberVariableGenes = config['tools']['analyse_citeseq']['numberVariableGenes'],
        outdir = 'results/citeseq_analysis/{sample}/' 
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:config['computingResources']['mediumRequirements']['threads']
    shell:
        "Rscript workflow/scripts/analyse_citeseq.R " +
        "--RDS {input.RDS}" +
        "--cellrangerADT {params.ADTFolder}" +
        "--h5 {input.h5}" +
        "--colorConfig {params.colorConfig}" +
        "--lookup {params.lookup}" +
        "--threshold {input.thresholds}" +
        "--threads {threads}" +
        "--sampleName {wildcards.sample}" +
        "--number_variable_genes {params.numberVariableGenes} " +
        "--output_directory {params.outdir} && touch {output.completeFile}"
