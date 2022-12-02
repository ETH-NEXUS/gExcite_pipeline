


#Rule to Analyse ADT Data in combination with GEXdata
rule plot_combined_ridgeplot:
    input:
        RDS = 'results/atypical_removed/{sample}.{demultiplexed}.genes_cells_filtered.corrected.atypical_removed.RDS' ,
        CellrangerADT = 'results/cellranger_adt/{sample}.{demultiplexed}.matrix.mtx',
    output:
        completeFile = 'results/cohort_combined_analysis/{sample}.plot_combined_ridgeplot.complete.txt'
    conda:
        "../envs/plot_combined_ridgeplot.yaml"
    params:
        ADTFolder = gatherCellrangerADTFolder() #'results/cellranger_adt/{sample}/outs/',
        outdir = 'results/cohort_combined_analysis/'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem_mb'],
        time_min = config['computingResources']['mediumRequirements']['runtime']
    threads:config['computingResources']['mediumRequirements']['threads']
    shell:
        "Rscript ../scripts/plot_citeseq_combined_ridgeplots.R " +
        "--atypicalRemoved {input.RDS}" +
        "--analysisADT {params.ADTFolder}" +
        "--sampleNames {wildcards.sample}" +
        "--outfolder {params.outdir} && touch {output.completeFile}"
