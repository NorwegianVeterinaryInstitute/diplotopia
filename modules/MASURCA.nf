process MASURCA {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //conda (params.enable_conda ? 'bioconda::chewbbaca=3.3.1' : null)
    
    container 'quay.io/biocontainers/masurca:4.1.1--pl5321hb5bd705_0'

    label 'process_high_memory'

    tag "$assemblyID" 
    
    input:
    tuple val(assemblyID), path(R1), path(R2), path(long_read)

    output: 
    tuple val(assemblyID), path("*"), emit: masurca_ch

    script: 
    """
    masurca -t 4 -i ${R1},${R2} -r ${long_read} 2>&1 | tee ${assemblyID}_masurca.log
    """
}