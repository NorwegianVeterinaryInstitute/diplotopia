// variant calling
process FREEBAYES_REGION {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //container 'quay.io/biocontainers/freebayes:1.3.7--h6a68c12_1' does not work 
    container 'quay.io/biocontainers/snippy:4.6.0--0'
    


    input:
    tuple path(ref), path(faidx)
    
    output:
    path("ref_regions.txt"), emit: regions
    path("*.version")

    script: 
    """
    freebayes --version > freebayes.version
    fasta_generate_regions.py $faidx $params.chunk_size > ref_regions.txt
    """
}

process FREEBAYES_CALL {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //container 'quay.io/biocontainers/freebayes:1.3.7--h6a68c12_1'
    container 'quay.io/biocontainers/snippy:4.6.0--0'

    label 'process_high_memory_cpu_time'
    tag "$ID" 

    input:
    tuple val(ID), val(ploidy), path(bam), path(bam_bai), path(ref), path(faidx), path(regions)
    
    output:
    tuple val(ID), path("*_unfiltered.vcf"), emit: vcf
    path("*")

    script: 
    
    // only runs if vcf_ref is provided (must be from same ref)
    if (params.vcf_ref != null) {
        """
        freebayes-parallel $regions $task.cpus -@ $params.vcf_ref --ploidy $ploidy \\
        --gvcf $params.freebayes_args -f $ref $bam > ${ID}_unfiltered.vcf

        freebayes --version > freebayes.version
        """

        
    }
    // runs by default
    else {
        """     
        freebayes-parallel $regions $task.cpus  --ploidy $ploidy  \\
        --gvcf $params.freebayes_args -f $ref $bam > ${ID}_unfiltered.vcf

        freebayes --version > freebayes.version
        """        

    }
}


process FREEBAYES_CALL_POP {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //container 'quay.io/biocontainers/freebayes:1.3.7--h6a68c12_1'
    container 'quay.io/biocontainers/snippy:4.6.0--0'

    label 'process_high_memory_cpu_time'
    tag "$IDS" 

    input:
    val(ID)
    path(bams)
    path(bais)
    tuple path(ref), path(faidx)
    path(regions)
    
    output:
    tuple val(IDS), path("pop_unfiltered.vcf"), emit: vcf
    path("*")

    script: 
    IDS = ID.join("_")
    //ploidy.unique() but must be identical 
    
    // only runs if vcf_ref is provided (must be from same ref)
    if (params.vcf_ref != null) {
        """
        freebayes-parallel $regions $task.cpus -@ $params.vcf_ref --ploidy 2 \\
        --gvcf $params.freebayes_args -f $ref $bams > pop_unfiltered.vcf

        freebayes --version > freebayes.version
        echo "Freebayes run with IDS: ${ID}" >> freebayes.log
        """

        
    }
    // runs by default
    else {
        """     
        freebayes-parallel $regions $task.cpus  --ploidy 2  \\
        --gvcf $params.freebayes_args -f $ref $bams > pop_unfiltered.vcf

        freebayes --version > freebayes.version
        echo "Freebayes run with IDS: ${ID}" >> freebayes.log
        """        

    }
}

// --genotype-qualities
// hwen ? before filter or
// https://github.com/brwnj/freebayes-nf/blob/master/main.nf
// bcftools norm -c all -f $fasta --multiallelics - --threads ${task.cpus} \
        // --output ${params.project}.vcf.gz --output-type z \
       //  ${params.project}_dirty_sorted.vcf.gz
       //   tabix -p vcf ${params.project}.vcf.gz

