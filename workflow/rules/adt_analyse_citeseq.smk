#Rule to create an initial threshold file
#Input: feature file from cellranger adt run
#Ouput: Threshold table containing the features in the first column and the sample name as a header for the second
# In the second column 0 are filled in. These can be changed later on for manual thresholding
rule create_initial_threshold_file:
    input:
        CellrangerADT = 'results/cellranger_adt/{sample}.features.tsv',
    output:
        thresholds = 'results/citeseq_analysis/{sample}.thresholds.tsv',
    resources:
        mem_mb = config['computingResources']['mem_mb']['low'],
        time_min = config['computingResources']['runtime']['low']
    log:
        "logs/create_initial_threshold_file/{sample}.log"
    benchmark:
        'results/citeseq_analysis/benchmark/{sample}.create_initial_threshold_file.benchmark'
    threads:config['computingResources']['threads']['low']
    shell:
        "echo -e 'Antibody\t{wildcards.sample}' > {output.thresholds} ; " +
	"awk ' {{ print $2 }}' {input.CellrangerADT} | sed 's/$/\t0/' >> {output.thresholds}"
        

#Rule to Analyse ADT Data in combination with GEXdata
#Input: RDS and h5 file from single sample GEX run, threshold file as generated by rule create_initial_threshold_file, Cellranger adt file as generated from rule Rscript_demultiplex_count_matrix
#Output: One folder under /results/citeseq_analysis/{sample} containing lots of diagnostic plots and an RDS file containing all the data.
rule Rscript_analyse_citeseq:
    input:
        RDS = 'results/atypical_removed/{sample}.atypical_removed.RDS' ,
        CellrangerADT = 'results/cellranger_adt/{sample}.matrix.mtx',
        thresholds = 'results/citeseq_analysis/{sample}.thresholds.tsv',
        h5 = 'results/counts_corrected/{sample}.corrected.variable_genes.h5',
    output:
        RDS = 'results/citeseq_analysis/{sample}/{sample}.GEX_cellrangerADT_SCE.RDS',
    conda:
        "../envs/adt_analyse_citeseq.yaml"
    params:
        ADTFolder = 'results/cellranger_adt/{sample}/outs/filtered_feature_bc_matrix',
        colorConfig = config['scampi']['resources']['colour_config'] ,
        lookup = config['resources']['adt_lookup'],
        numberVariableGenes = config['tools']['analyse_citeseq']['numberVariableGenes'],
        outdir = 'results/citeseq_analysis/{sample}/' 
    log:
        "logs/Rscript_analyse_citeseq/{sample}.log"
    benchmark:
        'results/citeseq_analysis/benchmark/{sample}.analyse_citeseq.benchmark'
    resources:
        mem_mb = config['computingResources']['mem_mb']['medium'],
        time_min = config['computingResources']['runtime']['medium']
    threads:config['computingResources']['threads']['medium']
    shell:
        "Rscript workflow/scripts/analyse_citeseq.R " +
        "--RDS {input.RDS} " +
        "--cellrangerADT {params.ADTFolder} " +
        "--h5 {input.h5} " +
        "--colorConfig {params.colorConfig} " +
        "--lookup {params.lookup} " +
        "--threshold {input.thresholds} " +
        "--threads {threads} " +
        "--sampleName {wildcards.sample} " +
        "--number_variable_genes {params.numberVariableGenes} " +
        "--output {params.outdir} &> {log}"
