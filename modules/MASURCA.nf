process MASURCA {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'quay.io/biocontainers/masurca:4.1.1--pl5321hb5bd705_0'

    label 'process_high_memory'

    tag "$meta.id" 
    
    input:
    tuple val(meta), path(R1), path(R2)
    path(long_reads), optional : true

    output: 
    tuple val(meta), path("${meta.id}/*"), emit: masurca_assembly
    // Check path how to give path("CA/") - need assembly and need all results in output ? 
    path "versions.yml", emit : versions

    script:
    def short_reads_param = "-i ${R1},${R2}"
    def long_reads_param = long_reads ? "-r ${long_reads}" : ""
    def log_file = "${meta.id}_masurca.log"

    """
    masurca -v > versions.yml

    masurca -t $task.cpus \\
        ${short_reads_param} \\
        ${long_reads_param} \\
        2>&1 | tee ${log_file} 
    """

}

// For masurca with config
// process MASURCA {
//     //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
//     //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
//     container 'quay.io/biocontainers/masurca:4.1.1--pl5321hb5bd705_0'

//     label 'process_high_memory'

//     tag "$assemblyID" 
    
//     input:
//     tuple val(assemblyID), path(R1), path(R2), path(long_read)

//     output: 
//     tuple val(assemblyID), path("*"), emit: masurca_ch

//     script: 
//     """
//     masurca -t 4 -i ${R1},${R2} -r ${long_read} 2>&1 | tee ${assemblyID}_masurca.log
//     """
// }