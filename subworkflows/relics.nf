    IDS = ID.join("_") 

        vcfgzs = 
        SAMTOOLS_INDEX_VAR.out.bam
        .map { (id) =  [it[0]] }
        .combine(INDEX_VT_NORMALIZE.out.vcfgz, by: 0)

    vcfgzs_ids = vcfgzs 
        .map { (id) =  [it[0]] }
        .collect()

    vcfgzs_vcfgz = vcfgzs
        .map { (vcfgz)  =  [it[1]] }
        .collect()

    vcfgzs_index = vcfgzs
        .map { (index)  =  [it[2]] }
        .collect()
    
    // channels ID by type - if we want 
    /* 

    id_hybrid_ch = 
        input_ch
        .filter {it[5] == "hybrid"}
        .map { (ID, assembly, type) =  [ it[0], it[1], it[5] ]}


    id_short_ch = 
        input_ch
        .filter {it[5] == "short"}
        .map { (ID, assembly, type) =  [ it[0], it[1], it[5] ]}

    id_long_ch =
        input_ch
        .filter {it[5] == "long"}
        .map { (ID, assembly, type) =  [ it[0], it[1], it[5] ]}

    */

/* Non paralell - its just crashing
process FREEBAYES_CALL {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    container 'quay.io/biocontainers/freebayes:1.3.7--h6a68c12_1'

    label 'process_high_memory_time'
    tag "$ID" 

    input:
    tuple val(ID), val(ploidy), path(bam), path(bam_bai)
    tuple path(ref), path(faidx)
    
    output:
    tuple val(ID), path("*_unfiltered.vcf"), emit: vcf
    path("*.version")

    script: 
    
    // only runs if vcf_ref is provided (must be from same ref)
    if (params.vcf_ref != null) {
        """    
        freebayes -@ $params.vcf_ref --ploidy $ploidy \\
        --gvcf $params.freebayes_args -f $ref $bam > ${ID}_unfiltered.vcf

        freebayes --version > freebayes.version
        """

        
    }
    // runs by default
    else {
        """     
        freebayes --ploidy $ploidy  \\
        --gvcf $params.freebayes_args -f $ref $bam > ${ID}_unfiltered.vcf

        freebayes --version > freebayes.version
        """        

    }
}


process FREEBAYES_CALL_POP {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    container 'quay.io/biocontainers/freebayes:1.3.7--h6a68c12_1'

    label 'process_high_memory_time'
    tag "$IDS" 

    input:
    val(ID), 
    path(bams)
    path(bais)
    tuple path(ref), path(faidx)
    
    output:
    tuple val(IDS), path("${IDS}_unfiltered.vcf"), emit: vcf
    path("*.version")

    script: 
    IDS = ID.join("_")
    
    // only runs if vcf_ref is provided (must be from same ref)
    if (params.vcf_ref != null) {
        """    
        freebayes -@ $params.vcf_ref --ploidy 2 \\
        --gvcf $params.freebayes_args -f $ref $bams > ${IDS}_unfiltered.vcf

        freebayes --version > freebayes.version
        """

        
    }
    // runs by default
    else {
        """     
        freebayes --ploidy 2 \\
        --gvcf $params.freebayes_args -f $ref $bams > ${IDS}_unfiltered.vcf
        
        freebayes --version > freebayes.version
        """        

    }
}

process FREEBAYES_CALL_POP2 {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    //container 'quay.io/biocontainers/freebayes:1.3.7--h6a68c12_2' is bugging
    //container 'quay.io/biocontainers/freebayes:1.3.6--hb0f3ef8_7'
    container 'quay.io/biocontainers/freebayes:1.3.7--h6a68c12_1'
    // update to 1.3.8 when ready 

    //label 'process_high_memory_cpu_time'
    //tag "$IDS" 

    input:
    // ID and bam are collections for all samples
    val(ID)
    // merged bam 
    tuple path(bam), path(bam_bai)
    tuple path(ref), path(faidx)

    output: 
    path("unfiltered.vcf"), emit: vcf
    path("ref_regions.txt")
    path("freebayes.log")
    path("*.version")

    script: 
    IDS = ID.join("_")
    
    if (params.vcf_ref != null) {
        // only runs if vcf_ref is provided (must be from same ref)

        """
        # common 
        freebayes --version > freebayes.version
        fasta_generate_regions.py $faidx $params.chunk_size > ref_regions.txt
        echo "IDS:\n $IDS" > freebayes.log

        # was weird variants 
        #freebayes-parallel ref_regions.txt $task.cpus -@ $params.vcf_ref --ploidy 2 \\
        #--gvcf $params.freebayes_args -f $ref $bam > unfiltered.vcf

        freebayes-parallel ref_regions.txt $task.cpus -@ $params.vcf_ref \\
        $params.freebayes_args -f $ref $bam > unfiltered.vcf

        
        echo "freebayes run vcf reference $params.vcf_ref" >> freebayes.log
        """

        
    }
    // --gvcf coverage information  when variant not called (should be then repeats mainly)
    // TODO see how I should use the ploidy option here !  -p $ploidy - if not equal ? ok as diploids though
    // runs by default
    else {
        """
        # common 
        freebayes --version > freebayes.version
        fasta_generate_regions.py $faidx $params.chunk_size > ref_regions.txt
        echo "IDS:\n $IDS" > freebayes.log

        freebayes-parallel ref_regions.txt $task.cpus  --ploidy 2  \\
        --gvcf $params.freebayes_args -f $ref $bam > unfiltered.vcf

        #freebayes -f $ref --ploidy 2 --gvcf $params.freebayes_args $bam  > unfiltered.vcf
        echo "freebayes run without reference vcf" >> freebayes.log
        """        

    }
}

// freebayes -f ref.fa -@ in.vcf.gz aln.bam >var.vcf (does it need to be compressed ? )
// TODO : some filters for guetting more efficient - but need first some experience I think
// TODO --min-supporting-mapping-quality
// TODO --min-supporting-base-quality
// --read-mismatch-limit
// TODO - can use that for paraelisation ? https://github.com/brwnj/freebayes-nf/blob/master/main.nf





*/

process FREEBAYES_NORMALIZE  {
    //conda "/cluster/projects/nn9305k/src/miniconda/envs/"
    container 'quay.io/biocontainers/freebayes:1.3.7--h6a68c12_2'
    // update to 1.3.8 when ready 

    label 'process_short'
    tag "$ID" 

    input:
    path(filteredvcf)
    

    output: 
    path("normalized.vcf"), emit: vcf_ch

    script: 
    threads = task.cpus *2 
    memperthread = "${( task.memory.toGiga() / threads ).trunc(0)}G"

    """
    freebayes $filteredvcf | vcfallelicprimitives -kg >calls.vcf
    vcffilter -f "QUAL > ${params.freebayes_qual}" > filtered.vcf
    """
}