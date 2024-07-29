// Mapping of nanopore reads to its own assembly - Indexing assembly - Indexing bam
process DOTPLOT {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/dgenies"
    // container 'quay.io/biocontainers/dgenies:1.5.0--pyhdfd78af_1'
    //dgenies is only interactive



    
    label ''
    
    input:
    path(json_short_summary)
    
    output: 
    path("busco_report*")
    
    script:
    """
    # ? busco plots also ? 
    quarto render busco_summary.qmd --to html --output_file busco_report.html --no-cache --dir . --save_dir .
    """
}