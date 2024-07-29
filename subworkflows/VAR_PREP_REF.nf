// preparation of reference for variant calling
// low complexity region masking, high coverage region masking ...

include { VAR_PREP_BASE as VAR_PREP_BASE_REF } from '../subworkflows/VAR_PREP_BASE.nf'
// maybe change VAR_PREP_BASE as VAR_PREP_BASE_REF as done for samples

include { BBMASK } from '../modules/BBTOOLS.nf'
include { SAMTOOLS_FAIDX } from '../modules/SAMTOOLS.nf'
include { SEQKIT_ASSEMBLY_STATS } from '../modules/SEQKIT.nf'


workflow VAR_PREP_REF {

    take:
    ref_input_ch // channel: tuple val(ID) path(assembly), path(R1), path(R2)
    
    main:  
    // note eventually create a bam file with reads to filter eg. regions of high coverage
    // if want to this for reference masking
    if (params.bbmask_maxcov != null) {
        VAR_PREP_BASE_REF(ref_input_ch)   

        /* 
        // !! NOT IMPLEMENTED YET 
        // Indexing ref
        SAMTOOLS_FAIDX(ref_input_ch.map { (ID, assembly) = [it[0], it[1]] }
        
        // keeping reads with high coverage for masking
        // chekc if ok still have the high coverage after removal duplicates !!!  
        SAMTOOLS_FILTER()
        // bam to sam and indexing 
        SAMTOOLS_TOSAM(SAMTOOLS_FILTER.out.map { (ID, bam) = [it[0], it[1]])
      
        
        // masking low complexity regions in reference 
        BBMASK(SAMTOOLS_TOSAM.out.sam.)
        */
        println("Not implemented yet")

    }
    else {
        dummy_sam_ch = Channel.fromPath("$baseDir/assets/dummy.sam", checkIfExists: true)
        // Masks low entropy regions in reference
        BBMASK(
            ref_input_ch.map { (ID, assembly) = [it[0], it[1]] },
            dummy_sam_ch)
    }
    // reindexing masked assembly
    SAMTOOLS_FAIDX(BBMASK.out.ref_ID)

    // getting stats for the masked reference
    SEQKIT_ASSEMBLY_STATS(BBMASK.out.ref_ID)

    emit:
    masked_ref = BBMASK.out.ref // channel : [ path(masked_ref) ] softmasked ref for variant calling 
    // should we change to 
    // ref_index = VAR_PREP_BASE_REF.ref_index  - or not better to use the samples one
    
}