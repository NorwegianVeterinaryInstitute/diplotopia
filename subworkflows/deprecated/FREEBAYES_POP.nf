// variant calling with freebayes subworflow
// could try to paralellize leater on : look at  https://github.com/brwnj/freebayes-nf
// and improve the writing as nf core style ...

include { SAMBAMBA_MARKDUP } from "../modules/SAMBAMBA.nf"
include { SAMTOOLS_COV; SAMTOOLS_FAIDX; SAMTOOLS_MERGE_VAR; SAMTOOLS_SORT_VAR } from "../modules/SAMTOOLS.nf"
include { FREEBAYES_CALL } from "../modules/FREEBAYES.nf"
// include { FREEBAYES_CALL; FREEBAYES_FILTER; FREEBAYES_NORMALIZE } from "../modules/FREEBAYES.nf"


workflow FREEBAYES_POP {
    take:
    bam // channel: [val(ID), val(ploidy), path(bam)]
    bam_bai // channel: [path(bam_bai)]
    fasta // channel: /path/to/ref.fasta  ref_ch

    
    main:

    //mark duplicates (or remove duplicates)
    if (params.markdup == "mark" | params.markdup == "remove" ) {
        SAMBAMBA_MARKDUP(bam)
        new_bam = SAMBAMBA_MARKDUP.out.bam
    } else if (params.markdup == null) {
        new_bam = bam
    } else {
        exit 1, "Please review the parameters for marking or removing duplicates reads;\
        must be either 'mark' or 'remove' or 'null'."
    }

    // getting coverage stat of the sampled (ID) mapped on the reference
    SAMTOOLS_COV(new_bam) 
    //SAMTOOLS_COV.out.coverage_stats // channel: [val(ID), path(coverage_data)]

    // create index again - as it does not follow
    SAMTOOLS_FAIDX(fasta)

    // MODULE call variants - population of samples 
    // TODO Right now does not make use ploidy -> either remove from params and assume is 2 
    // TODO OR find solution to use it better 
    

    POPULATION VARIANT CALING - Not working as it should right now 
    // Collect all IDs and all bam files

    bams = new_bam
        .map { bam  =  it[2] }
        .collect()

    ids = new_bam
        .map { id =  it[0] }
        .collect()

    //SAMTOOLS_FAIDX.out.faidx_ch.view() // channel: [path(assembly), path(fai) ]


    // --- Calling per population - all isolates at once 
    // merging all the bam files into one - seems otherwise bug in parallel
    SAMTOOLS_MERGE_VAR (bams, bam_bai.collect())
    SAMTOOLS_SORT_VAR(SAMTOOLS_MERGE_VAR.out.bam) 

    // This will be one run only so ok to implement as is
    FREEBAYES_CALL_POP(ids, SAMTOOLS_SORT_VAR.out.bam, SAMTOOLS_FAIDX.out.faidx_ch) 
    // MODULE filter variants
    //FREEBAYES_FILTER()
    // MODULE normalize variants
    //FREEBAYES_NORMALIZE()

   
    
    //emit:
    //call_vcf = FREEBAYES_CALL.out.vcf_ch // channel: [val(ID), path(vcf)]
    //filtered_vcf = FREEBAYES_FILTER.out.vcf_ch // channel: [val(ID), path(vcf)]
    //normalized_vcf = FREEBAYES_NORMALIZE.out.vcf_ch // channel: [val(ID), path(vcf)]
    // to do after nf core way 
    //freebayes_version == FREEBAYES_CALL.version // path *.version 
    // stats =     
    
}


