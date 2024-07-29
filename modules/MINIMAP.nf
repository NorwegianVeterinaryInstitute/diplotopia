// Mapping of nanopore reads to its own assembly - Indexing assembly - Indexing bam
process MINIMAP {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/minimap_2.26"
    container 'quay.io/biocontainers/minimap2:2.26--he4a0461_2'

    tag "$ID"
    
    label 'process_high'
    
    input:
    tuple val(ID), path(assembly), path(longreads)
    
    output: 
    tuple val(ID), path("${ID}_long.sam"), emit: longsam_ch
    
    script:
    threads = task.cpus *2
    """
    minimap2 -d ${assembly}.mmi $assembly 
    minimap2 -t $threads -ax map-ont ${assembly}.mmi $longreads -o ${ID}_long.sam
    minimap2 --v > minimap2.version
    """
}

// see if find a way to merge with above - but output must be in - define out separately ? 
process HAPLO_MINIMAP {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/minimap_2.26"
    container 'quay.io/biocontainers/minimap2:2.26--he4a0461_2'

    tag "$ID"
    
    label 'process_high'
    
    input:
    tuple val(ID), path(assembly), path(longreads)
    
    output: 
    tuple val(ID), path(assembly), path("${ID}_long.sam"), emit: longsam_ch
    
    script:
    threads = task.cpus *2

    """
    minimap2 -d ${assembly}.mmi $assembly 
    minimap2 -t $threads -ax map-ont ${assembly}.mmi $longreads --secondary=no -o ${ID}_long.sam
    minimap2 --v > minimap2.version
    """
}