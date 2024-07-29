// Kraken2 on assemblies
process KRAKEN {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/kraken2"
    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0'
    
    //label 'process_high'
    tag "$ID"
    
    input:
    tuple val(ID), path(assembly)
    
    output: 
    path("*")


    script:
    """
    kraken2 -db $params.krakenDB \\
    --threads $task.cpus \\
    --use-names \\
    --output  ${ID}.out \\
    --report  ${ID}.report \\
    $assembly

    kraken2 -v > kraken2.version
    """
}

