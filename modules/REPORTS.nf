// Mapping of nanopore reads to its own assembly - Indexing assembly - Indexing bam
process BUSCO_REPORT {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/R"
    container 'evezeyl/compassr:latest'

    tag "Busco_summary"
    
    label 'process_medium_memory'
    
    input:
    path(json_short_summary)
    
    output: 
    path("*")
    
    script:
    """
    cp $baseDir/bin/busco_summary.qmd .

    quarto render busco_summary.qmd \\
    -P res_dir:"." \\
    -P save_dir:"." \\
    -P lineage_dataset:$params.lineage_dataset 2>&1 | tee BUSCO_SUMMARY.log

    # copy the .command.sh to output if need adjustement
    cp .command.sh busco_summary.command.sh
    """
}