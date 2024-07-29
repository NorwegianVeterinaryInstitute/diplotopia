//concatenate reads when run in several files 

process CONCAT_READS { 

    input:
    tuple val(assemblyID), path(samfile), path(assembly)

    output: 
    tuple val(ID), path("*.bam"), emit: concat_reads_ch
    // could be better solution
    // https://stackoverflow.com/questions/74039553/nextflow-rename-barcodes-and-concatenate-reads-within-barcodes

    script: 
    """
    """