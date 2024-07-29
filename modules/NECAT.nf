// preprocessing module & 
// correction module: progressive strategy to correct nanopore reads (2 steps)
// output nbX correction raw reads ./$PROJECT/1-consensus/cns_iter${NUM_ITER}/cns.fasta
// output nbX correction raw reads ./$PROJECT/1-consensus/cns_final.fasta

process NECAT_GLOBAL {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    container 'evezeyl/necat:latest'

    tag "$assemblyID"

    label 'process_high_memory'

    input:
    // this should link the readlist to the path 
    tuple val(assemblyID), path(necatconf), path(readlist)

    output:
    tuple val(assemblyID), path(necatconf), path(readlist), path("**"), emit: necat_ch
   

    script: 
    """
    necat.pl correct ${necatconf} 2>&1 | tee ${assemblyID}_necat_correct.log
    necat.pl assemble ${necatconf} 2>&1 | tee ${assemblyID}_necat_assemble.log
    necat.pl bridge ${necatconf} 2>&1 | tee ${assemblyID}_necat_bridge.log
    """
}

// --------------------------------------------------------
// Those should work but it seems some results are lacking
// --------------------------------------------------------

// trimming module (I guess &)
// assembly module : builds string graph to assemble genome (2 steps)

process NECAT_ASSEMBLE {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'evezeyl/necat:latest'

    label 'process_medium_memory'

    tag "$assemblyID" 
    
    input:
    tuple val(assemblyID), path(necatconf), path("*")

    output: 
    // to modify to get the NECAT project config with assemblyID and then path
    tuple val(assemblyID), path(necatconf), path("*"), emit: necat_assemble_ch

    script: 
    """
    necat.pl assemble ${necatconf} 2>&1 | tee ${assemblyID}_necat_assemble.log
    """
}

// bridging & polishing if option set
process NECAT_BRIDGE {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'evezeyl/necat:latest'

    label 'process_medium_memory'

    tag "$assemblyID" 
    
    input:
    tuple val(assemblyID), path(necatconf), path("*")

    output: 
    // to modify to get the NECAT project config with assemblyID and then path
    tuple val(assemblyID), path(necatconf), path("*"), emit: necat_bridge_ch

    script: 
    """
    necat.pl bridge ${necatconf} 2>&1 | tee ${assemblyID}_necat_bridge.log
    """
}

// preprocessing module & 
// correction module: progressive strategy to correct nanopore reads (2 steps)
// output nbX correction raw reads ./$PROJECT/1-consensus/cns_iter${NUM_ITER}/cns.fasta
// output nbX correction raw reads ./$PROJECT/1-consensus/cns_final.fasta
process NECAT_CORRECT {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    container 'evezeyl/necat:latest'

    tag "$assemblyID"

    label 'process_medium_memory'

    
    input:
    tuple val(assemblyID), path(necatconf)

    output:
    // to modify to get the NECAT project config with assemblyID and then path
    tuple val(assemblyID), path(necatconf), path("*"), emit: necat_correct_ch
   

    script: 
    """
    necat.pl correct ${necatconf} 2>&1 | tee ${assemblyID}_necat_correct.log
    """
}