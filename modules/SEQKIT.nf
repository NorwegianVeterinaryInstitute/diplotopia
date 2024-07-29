//SPLIT each assembly into contigs - prepare for BLASTING
process SEQKIT_TOCONTIG { 

   	//conda (params.enable_conda ? 'bioconda::bwa=0.7.8' : null)
	container 'quay.io/biocontainers/seqkit:2.8.1--h9ee0642_0'


    tag "${ID}_${contigID}"

    input:
    tuple val(ID), path(assembly)

    output: 
    tuple val(ID), path("*.part_*"), emit: contigs_ch
    path("*")

    script:
    """
    # prepare a list of correspondance in case
    grep ">" $assembly > ${ID}_contigs_ids.txt 

    seqkit split $assembly -i -t dna \\
    -O .  \\
    --threads $task.cpus

    seqkit version > seqkit.version
    """
}

//SPLIT each assembly into contigs - prepare for BLASTING
process SEQKIT_ASSEMBLY_STATS { 

   	//conda (params.enable_conda ? 'bioconda::bwa=0.7.8' : null)
	container 'quay.io/biocontainers/seqkit:2.8.1--h9ee0642_0'


    tag "${ID}"
    
    input:
    tuple val(ID), path(assembly)

    output: 
    tuple val(ID), path("*assembly_stats.tsv"), emit: seqkit_assembly_stats

    script:
    """
    # put the gap letter to N to compute masking
    seqkit stats -a -G "N" $assembly > ${ID}_seqkit_assembly_stats.tsv
    seqkit version > seqkit.version
    """
}


