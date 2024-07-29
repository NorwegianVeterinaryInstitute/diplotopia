// used in haplopurge track
process HAPLOPURGE_PURGE_HIST {

    conda "/cluster/projects/nn9305k/src/miniconda/envs/purge_haplotigs"
    container 'quay.io/biocontainers/purge_haplotigs:1.1.3--hdfd78af_0'

    tag "$ID"
    
    label 'process_short'


    input:
    tuple val(ID), path(assembly), path(bam)
    
    
    output: 
    tuple val(ID), path(assembly), path("*.gencov"), emit: haplohist_ch
    path("*")

    script:
    // need to path the index files is prefered? 
    """
    purge_haplotigs  hist  -b $bam  -g $assembly -d $params.purgehaplotigs_maxdepth -t $task.cpus
    mv tmp_purge_haplotigs ${ID}_tmp_purge_haplotigs    
    """
}
// used in variant track - can we merge those ? 
process PURGE_HIST {

    conda "/cluster/projects/nn9305k/src/miniconda/envs/purge_haplotigs"
    container 'quay.io/biocontainers/purge_haplotigs:1.1.3--hdfd78af_0'

    tag "$ID"
    
    label 'process_short'


    input:
    // ! to modify was origin in haplopurge
    //tuple val(ID), path(assembly), path(bam)
    tuple val(ID), path(assembly), path(index), path(bam), path(bai)
    
    output: 
    tuple val(ID), path(assembly), path("*.gencov"), emit: haplohist_ch
    path("*")

    script:
    // need to path the index files is prefered? 
    """
    purge_haplotigs  hist  -b $bam  -g $assembly -d $params.purgehaplotigs_maxdepth -t $task.cpus
    mv tmp_purge_haplotigs ${ID}_tmp_purge_haplotigs    
    """
}

// those params need to move into the csv input file at one point if want individual parameters
// for multiple isolates 

process PURGE_COVERAGE {

    conda "/cluster/projects/nn9305k/src/miniconda/envs/purge_haplotigs"
    container 'quay.io/biocontainers/purge_haplotigs:1.1.3--hdfd78af_0'

    tag "$ID"
    label 'process_short'


    input:
    tuple val(ID), path(assembly), path(gencov), val(haplo_low), val(haplo_mid), val(haplo_high)
    
    output: 
    tuple val(ID), path(assembly), path("*_coverage_stats.csv"), emit: haplopurged_ch
    path("*")

    script:
    """
    purge_haplotigs  cov  -i $gencov  -l $haplo_low  -m $haplo_mid  -h $haplo_high  \\
    -o ${ID}_coverage_stats.csv -j $params.haplo_junk  -s $params.haplo_suspect 
    """
}

// when do for all : need mv tmp_purge_haplotigs ${ID}_tmp_purge_haplotigs - so we can keep all directories

process PURGE_HAPLOTIGS {

    conda "/cluster/projects/nn9305k/src/miniconda/envs/purge_haplotigs"
    container 'quay.io/biocontainers/purge_haplotigs:1.1.3--hdfd78af_0'

    tag "$ID"
    label 'process_high_memory_time'

    input:
    tuple val(ID), path(assembly), path(coverage_stats), path(bam), path(index)
    
    output: 
    tuple val(ID), path("*_haplocurated.fasta"), emit: haplopurged_ch
    path("*")
    // channel emit for chaining eventually 

    script:
    threads = task.cpus * 2 
    """
    purge_haplotigs purge -t $threads -o ${ID}_haplocurated -d -b $bam \\
    -a $params.haplo_align_cov -m $params.haplo_max_match \\
    -g $assembly  -c $coverage_stats

    mv tmp_purge_haplotigs ${ID}_tmp_purge_haplotigs
    """
}
