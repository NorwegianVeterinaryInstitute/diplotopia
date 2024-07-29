// Base for variant calling with freebayes subworflow

// for both 
include { FREEBAYES_REGION } from "../modules/FREEBAYES.nf"
include { BCFTOOLS_FILTER } from "../modules/BCFTOOLS.nf"
include { VT_NORMALIZE } from "../modules/VT.nf"

include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FREEBAYES_CALL;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_BCFTOOLS_FILTER;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_VT_NORMALIZE } from "../modules/BCFTOOLS.nf"


// for samples
include { FREEBAYES_CALL } from "../modules/FREEBAYES.nf"

// for population 
include { FREEBAYES_CALL_POP } from "../modules/FREEBAYES.nf"

// for combining from individual calls 
include { BCFTOOLS_MERGE } from "../modules/BCFTOOLS.nf"



workflow CALL_SAMPLES {
    take:
    bam_bai // bam_bai -> if not indexed need SAMTOOLS_INDEX_VAR 
    ploidy // channel with ploidy val(ID), val(ploidy) 
    ref_index // path (ref) path (index)
    
    
    main:
    // Calling on individual samples    
    // Region split of reference - for parallelisation  
    FREEBAYES_REGION(ref_index)
    
    // Variant calling per sample 
    //  tuple val(ID), val(ploidy), path(bam), path(bam_bai), path(ref), path(faidx), path(regions)
    FREEBAYES_CALL(
        ploidy
        .combine(bam_bai, by : 0)
        .combine(ref_index)  
        .combine(FREEBAYES_REGION.out.regions)
        )

    
    // --------------------------------
    // Filter variants 
    BCFTOOLS_FILTER(FREEBAYES_CALL.out.vcf)
    // Normalize variants 
    VT_NORMALIZE(BCFTOOLS_FILTER.out.vcf
                 .combine(ref_index.map{ (ref) = [it[0]] })
                ) 


    // Indexing the different types of vairants (for IGV visualisation)
    BCFTOOLS_INDEX_FREEBAYES_CALL(FREEBAYES_CALL.out.vcf)
    BCFTOOLS_INDEX_BCFTOOLS_FILTER(BCFTOOLS_FILTER.out.vcf)
    BCFTOOLS_INDEX_VT_NORMALIZE(VT_NORMALIZE.out.vcf) 
    
    
    
    emit:
    // vcf_unfiltered = FREEBAYES_CALL.out.vcf  // tuple val(ID), path("*_unfiltered.vcf")
    // vcf_filtered =BCFTOOLS_FILTER_SAMPLES.out.vcf  // tuple val(ID), path("*_filtered.vcf")
    norm_vcfgz_index = BCFTOOLS_INDEX_VT_NORMALIZE.out.vcfgz_index  // tuple val(ID), path("*_normalized.vcf.gz"), path(*tbi) - normalized variants



}

workflow CALL_POP {

    take:
    pop_id // collection of all ids  
    pop_bam // collection of all bams 
    pop_bai // collection of all bai 
    ref_index // path reference,  path index

    // only possible for diploids now default is 

    main:
    // Region split of reference - for parallelisation  (make it also here - if want to do only population)
    FREEBAYES_REGION(ref_index)
    // variant call as a population 
    // TODO check be sure that ploidy 2 is what to put here ! on the merge and not the total ploidy
    FREEBAYES_CALL_POP(pop_id, pop_bam, pop_bai, ref_index, FREEBAYES_REGION.out.regions)

    // --------------------------------
    // Filter variants 
    BCFTOOLS_FILTER(FREEBAYES_CALL_POP.out.vcf)

    // Normalize variants 
    VT_NORMALIZE(BCFTOOLS_FILTER.out.vcf
                 .combine(ref_index.map{ (ref) = [it[0]] })
                ) 

    // Indexing the different types of vairants (for IGV visualisation)
    BCFTOOLS_INDEX_FREEBAYES_CALL(FREEBAYES_CALL_POP.out.vcf)
    BCFTOOLS_INDEX_BCFTOOLS_FILTER(BCFTOOLS_FILTER.out.vcf)
    BCFTOOLS_INDEX_VT_NORMALIZE(VT_NORMALIZE.out.vcf) 
    
    emit:
    // vcf_unfiltered = FREEBAYES_CALL_POP.out.vcf  // tuple val(IDS), path("*_unfiltered.vcf")
    // vcf_filtered =BCFTOOLS_FILTER_POP.out.vcf  // tuple val(ID), path("*_filtered.vcf")
    pop_norm_vcfgz_index = BCFTOOLS_INDEX_VT_NORMALIZE.out.vcfgz_index  // tuple val(ID), path("*_normalized.vcf.gz"), path(*tbi) - normalized variants
    
}



// Combine individual variant calls - for crosscheck with population call 
// does not seems to work if uncompressed - so only do a the end 
// here would have better have its by index rearranged but did not manage now

workflow COMBINE_INDIV_CALLS {
    take:
    norm_vcfgz_index  // tuple val(ID), path("*_normalized.vcf.gz"), path(*tbi) - normalized variants  
    ref // path ref
    
    main:
    coll_id = norm_vcfgz_index
                .map { (id) =  [it[0]] }
                .collect()

    coll_vcfgz = norm_vcfgz_index
                .map { (vcfgz) =  [it[1]] }
                .collect()

    coll_tbi = norm_vcfgz_index
                .map { (tbi) =  [it[2]] }
                .collect()

    // merging individually called - to check with population 
    BCFTOOLS_MERGE(coll_id, coll_vcfgz, coll_tbi, ref)

    emit:
    combined_norm_vcfgz_index = BCFTOOLS_MERGE.out.vcfgz  // tuple val(IDS), path("*_merged.vcfgz"), path(_merged.vcfgz.tbi) - merged normalized variants
}