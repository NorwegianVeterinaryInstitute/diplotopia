// Preparing samples for variant calling with freebayes

// We need in term to be able to specify options ... as this 
/*
options.args = [ 
    "--db ${params.abricate_db}",
    "--minid ${params.minid}",
    "--mincov ${params.mincov}"
    ]
*/

include { VAR_PREP_BASE as VAR_PREP_BASE_SAMPLES } from '../subworkflows/VAR_PREP_BASE.nf'
include { SAMTOOLS_MERGE_BAMS; SAMTOOLS_COORDSORT_INDEX } from '../modules/SAMTOOLS.nf'

workflow VAR_PREP_SAMPLES {

    take: 
    illumina_input_ch // channel: [ID, R1, R2]
    masked_ref // channel : val(ID), path(ref) : softmasked ref for variant calling 
    


    main:
    // preparation for all samples individually - to the same reference
    // Reformulated from BASE - was problem to multiply the ref.

    VAR_PREP_BASE_SAMPLES(
        illumina_input_ch
        .combine(masked_ref)
        .map{(ID, assembly, R1, R2) = [it[0], it[3], it[1], it[2]] }
        )
        

    // -------------------------------------------------------------
    // Samples into population  
    // TODO make it optional ? 
    // TODO there is many duplicate files - file name colision on the output collision why ? 
    // check if we get the bam_bai channel
        
    pop_id = VAR_PREP_BASE_SAMPLES.out.bam_bai
        .map { (id) =  [it[0]] }
        .collect()

    pop_bam = VAR_PREP_BASE_SAMPLES.out.bam_bai
        .map { (bam)  =  [it[1]] }
        .collect()
    //pop_bam.view()

    pop_bai = VAR_PREP_BASE_SAMPLES.out.bam_bai
        .map { (bai)  =  [it[2]] }
        .collect()

    // merging all the bam files into one
    /// there is a name colision but it should not why 
    SAMTOOLS_MERGE_BAMS(pop_id, pop_bam, pop_bai)
    SAMTOOLS_COORDSORT_INDEX(SAMTOOLS_MERGE_BAMS.out.bam) 


    emit:
    bam_bai       = VAR_PREP_BASE_SAMPLES.out.bam_bai  // tuple: [ID, path(bam), path(bai)] 
    pop_id // collection of all ids  
    pop_bam // collection of all bams 
    pop_bai // collection of all bai 
    ref_index = VAR_PREP_BASE_SAMPLES.out.ref_index // path (ref) path (index) - just to follow 
   


}

