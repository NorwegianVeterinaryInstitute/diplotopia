process QUAST_REF {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/quast_5.2.0"
    container 'evezeyl/quast:latest'
    
    label 'process_high_memory_time'
    
    input:
    path(assemblies)
    
    output: 
    path("quast_ref/*")


    script:
    """
    quast --large \\
    -r $params.quast_ref \\
    --fragmented \\
    --features $params.quast_annot \\
    -o quast_ref \\
    --min-contig $params.min_contig \\
    --report-all-metrics --plots-format $params.plot_format \\
    --no-sv --threads $task.cpus $assemblies

    quast --version > quast.version
    """
}

// --gene-finding --eukaryote \\
// --rna-finding --conserved-genes-finding \\

// Use own assembly as reference

process QUAST_ASSEMBLY {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/quast_5.2.0"
    container 'evezeyl/quast:latest'
    
    label 'process_high_memory_time'

    input:
    path(assemblies)

    output: 
    path("assembly_ref/*")


    script:
    """
    quast --large \\
    -r $params.quast_ref_assembly \\
    --fragmented \\
    -o assembly_ref \\
    --min-contig $params.min_contig \\
    --report-all-metrics --plots-format $params.plot_format \\
    --no-sv --threads $task.cpus $assemblies
    
    quast --version > quast.version
    """
}

//  --gene-finding --eukaryote \\
//  --rna-finding --conserved-genes-finding \\
