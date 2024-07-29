// manipulating and re-indexing short reads 
// adding tag and process at rerun - choose which one to do
// tag "$ID" 
// label 'process_high'

// Adjust name
// adjust threads as SAMTOOLS_MERGE_VAR 
// threads = task.cpus * 2 - 1
// memperthread = "${( task.memory.toGiga() / threads ).trunc(0)}G"

process SAMTOOLS_INDEX_SHORT {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    input:
    tuple val(ID), path(shortbam)
 
    output: 
    tuple val(ID), path("${ID}_short.bam"), emit: shortindexbam_ch // to remove I think
    tuple val(ID), path("*.bam.bai"), emit: bai

    script: 
    """
    # reactualize index in case only using short reads
    samtools sort -o ${ID}_short.bam $shortbam
    # reactualize index
    samtools index ${ID}_short.bam

    samtools --version > samtools.version
    """
}


// Merging short and long reads that were mapped into one file
process SAMTOOLS_MERGE_LONG_SHORT {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0' 

    label 'process_high'
    
    tag "$ID" 

    input:
    tuple val(ID), path(longreadsbam), path(shortreadsbam)

    output: 
    tuple val(ID), path("${ID}.bam"), emit: mergedbam_ch

    script: 
    """   
    # merging short and long reads
    samtools merge merged.bam $longreadsbam $shortreadsbam
    
    # sorting by coordinate
    samtools sort -o ${ID}.bam merged.bam
    
    # reactualize index
    samtools index ${ID}.bam
    
    samtools --version > samtools.version
    """
}

// for several bams
// rewritten from SAMTOOLS_MERGE_VAR 
process SAMTOOLS_MERGE_BAMS {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0' 

    tag "$IDS" 
    label 'process_high_cpu_short'

    input:
    val(ID)
    path(bam)
    path(bai)
    
    output: 
    tuple val(IDS), path("*temp.bam"), emit: bam

    script: 
    IDS = ID.join("_")
    """
    samtools merge -@ $task.cpus -o {IDS}_merged_temp.bam $bam
    samtools --version > samtools.version
    """
}
// SAMTOOLS_SORT_VAR replaced by SAMTOOLS_CORDSORT_INDEX
// TODO adjust threads mem 

/* deprecated I think
process SAMTOOLS_COV {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0' 

    label 'process_high'
    
    tag "$ID" 

    input:
    tuple val(ID), val(ploidy), path(bam)
    
    output: 
    tuple val(ID), path("*data.tsv"), emit : coverage_stats
    path("*")

    script: 
    """ 
    samtools coverage --plot-depth --excl-flags UNMAP,SECONDARY,QCFAIL,DUP -o ${ID}_coverage_filtered_plot.tsv $bam 
    samtools coverage --excl-flags UNMAP,SECONDARY,QCFAIL,DUP -o ${ID}_coverage_filtered_data.tsv $bam 
    samtools coverage --plot-depth -o ${ID}_coverage_unfiltered_plot.tsv $bam 
    samtools coverage -o ${ID}_coverage_unfiltered_data.tsv $bam 

    samtools --version > SAMTOOLS_COV_samtools.version
    """
}
*/ 
// TODO : for now let the plots as is - though should do better looking plots maybe as summary 
// TODO : pending - as might change the software later on

process SAMTOOLS_FAIDX {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    fair true

    input:
    tuple val(ID), path(assembly)
 
    output: 
    tuple val(ID), path(assembly), path("*.fai"), emit: ref_index // modifiy for compass
    tuple path(assembly), path("*.fai"), emit: noid_assembly_fai

    script: 
    """
    samtools faidx $assembly    
    samtools --version > samtools.version
    """
}

// TODO ensure not deprecated ?
process SAMTOOLS_INDEX_VAR {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID" 

    input:
    tuple val(ID), val(ploidy), path(bam)
 
    output: 
    tuple val(ID), val(ploidy), path(bam), path("*.bam.bai"), emit: ploidy_bam_bai

    script: 
    """
    samtools --version > samtools.version
    samtools index $bam
    
    """
}

// VAR modules reqriten - modify in compass ?
process SAMTOOLS_TOBAM {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID"
    label 'process_high'

    input:
    tuple val(ID), path(sam)

    output: 
    tuple val(ID), path("*.bam"), emit: bam
    

    script: 
    """
    # transformation to bam
    samtools view -@ $task.cpus $sam -o ${ID}_map1.bam
    samtools --version > samtools.version
    """
}
// remove S option useless now 
process SAMTOOLS_TOBAM_LONG {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID"
    label 'process_high'

    input:
    tuple val(ID), path(longsam)

    output: 
    tuple val(ID), path("${ID}_long.bam"), emit: longbam_ch
    path("*.bai")

    script: 
    threads = task.cpus *2 
    memperthread = "${( task.memory.toGiga() / threads ).trunc(0)}G"

    """
    samtools --version > TOBAM_samtools.version

    samtools sort -m $memperthread -@ $threads -o sorted.sam $longsam
    samtools view -m $memperthread -@ $threads -bS sorted.sam -o ${ID}_long.bam
    samtools index ${ID}_long.bam
    samtools --version > samtools.version
    """
}

process HAPLO_TOBAM {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    
    label 'process_short'
    tag "$ID" 

    input:
    tuple val(ID), path(assembly), path(longsam)

    output: 
    tuple val(ID), path(assembly), path("*_long.bam"), emit: longbam_ch
    tuple val(ID), path("*_long.bam"), path("*_long.bam.bai"), emit: mapindex_ch 


    script: 
    threads = task.cpus *2 
    memperthread = "${( (task.memory.toGiga() -1) / threads ).trunc(0)}G"

    """
    samtools sort -m $memperthread -@ $threads -o sorted.sam $longsam
    samtools view -m $memperthread -@ $threads -b sorted.sam -o ${ID}_long.bam
    samtools index ${ID}_long.bam
    samtools --version > samtools.version
    """
}

process SAMTOOLS_NAMESORT {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID"
    label 'process_high'

    input:
    tuple val(ID), path(bam)

    output: 
    tuple val(ID), path("*_sorted.bam"), emit: bam

    script: 
    """
    # sorting by Name 
    samtools sort -@ $task.cpus -N -o ${ID}_sorted.bam $bam
    samtools --version > samtools.version
    """
}

process SAMTOOLS_FIXMATE {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID"
    label 'process_high'

    input:
    tuple val(ID), path(bam)

    output: 
    tuple val(ID), path("*_fixmate.bam"), emit: bam

    script: 
    """
    # fills mate coordinates and insert size - annotation map
    # -r remove secondary and unmapped reads
    # add ms score tags
    samtools fixmate -@ $task.cpus -r -m -c -O bam $bam ${ID}_fixmate.bam
    samtools --version > samtools.version
    """
}

process SAMTOOLS_COORDSORT_INDEX {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID"
    label 'process_high'

    input:
    tuple val(ID), path(bam)

    output: 
    tuple val(ID), path("*_coordsort.bam"), path("*bai"), emit: bam_bai

    script: 
    """
    # sorting alignment by coordinates - reactualize index
    samtools sort -@ $task.cpus -o ${ID}_coordsort.bam $bam
    samtools index -@ $task.cpus ${ID}_coordsort.bam
    samtools --version > samtools.version
    """
}

process SAMTOOLS_MARKDUP {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID"
    label 'process_high'

    input:
    tuple val(ID), path(bam), path(bai)

    output: 
    tuple val(ID), path("*_markdup.bam"), path("*.bai"), emit: bam_bai
    path("*_markdup.stats")
    path("*")

    script: 
    """
    # marks dupplicate alignments; mark supplementary reads as duplicates 
    # removes duplicate reads 
    # in: sorted fixed mate coord indexed bam 
    samtools markdup -@ $task.cpus -r -S \\
    --duplicate-count -m s \\
    --use-read-groups \\
    -s -f ${ID}_markdup.stats \\
    $bam ${ID}_markdup.bam
    samtools --version > samtools.version

    # need to reactualize index to be sure
    samtools index -@ $task.cpus ${ID}_markdup.bam
    
    samtools --version > samtools.version
    """
}
// https://www.htslib.org/doc/samtools-stats.html
// https://stackoverflow.com/questions/54908223/plot-bamtools-where-are-the-commands
// Might replace flagstat 
process SAMTOOLS_BAMSTATS {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID"
    label 'process_high'

    input:
    tuple val(ID), path(bam), path(bai)
    
    output:
    tuple val(ID), path("*.bamstats"), emit: bamstats
    
    script:
    """
    # gets bam statistics - bamfile must be sorted 
    samtools stats -@ $task.cpus $bam > ${ID}.bamstats 
    # seems cannot create plot here 
    samtools --version > samtools.version
    """

}

process SAMTOOLS_FLAGSTAT {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID"
    label 'process_high'

    input:
    tuple val(ID), path(bam), path(bai)
    
    output:
    tuple val(ID), path("*_flags.stats"), emit: flagstats
    
    script:
    """
    # gets flags statistics
    samtools flagstat -@ $task.cpus $bam -O tsv > ${ID}_flags.stats 
    samtools --version > samtools.version
    """

}

process SAMTOOLS_FILTER {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID"
    label 'process_high'

    input:
    // assembly is the ref used for mapping and its index
    tuple val(ID), path(bam), path(bam_bai), path(assembly), path(index)
    
    output:
    tuple val(ID), path("*_map2.bam"), path("*_map2.bam.bai"), emit: bam_bai
    tuple val(ID), path(assembly), path(index), path("*_map2.bam"), path("*_map2.bam.bai"), emit : bam_tohist
    path("*.unmapped")


    script:
    """
    # filtering by quality, unmapped, secondary, supplementary, dup
    samtools view -@ $task.cpus -b -h \\
    -q $params.min_mapq $params.filterout_tags \\
    --sanitize all -T $assembly -t $index \\
    -U ${ID}.unmapped -o ${ID}_map2.bam $bam

    # reindexing by coordinates
    samtools index -@ $task.cpus ${ID}_map2.bam
    samtools --version > samtools.version
    """

}

process SAMTOOLS_COVERAGE {
    conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost2022-bwa"
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    tag "$ID"
    label 'process_high'

    input:
    tuple val(ID), path(bam), path(bai)
    
    output:
    tuple val(ID), path("*.depth"), path("*.depth_all_pos"), emit: coverage
    path("*")
    
    script:
    """
    # compute depth at each position - incl. gaps
    samtools depth -J -a $bam -o ${ID}.depth_all_pos
    samtools depth -J $bam -o ${ID}.depth

    samtools coverage --plot-depth -o ${ID}_coverage.plot $bam
    samtools coverage -o ${ID}_coverage.data $bam

    samtools --version > samtools.version
    """

}


/*     

    # 


 
 */ 

/*
FLAGS 
0x0004	u	the query sequence itself is unmapped  = UNMAP 
0x0008	U	the mate is unmapped = MUMAMP
0x0100	s	the alignment is not primary (secondary) = SECONDARY
0x0400	d	the read is either a PCR or an optical duplicate = DUP 
0x0800	S	the alignment is supplementary = SUPPLEMENTARY

0x4000  
flags QCFAIL, 

*/

// Samtools stats https://www.biostars.org/p/9545575/
