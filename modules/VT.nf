process VT_NORMALIZE {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    // https://genome.sph.umich.edu/wiki/Vt#Introduction

    container 'quay.io/biocontainers/snippy:4.6.0--0'

    label 'process_high_memory_time2'
    tag "$ID"

    input:
    tuple val(ID), path(filtered), path(ref)

    output: 
    tuple val(ID), path("*normalized.vcf"), emit: vcf
    //path("*")

    script: 
    """
    vt normalize $filtered -r $ref -o ${ID}_normalized.vcf
    """

}