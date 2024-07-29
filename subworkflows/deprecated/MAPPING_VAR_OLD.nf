// Mapping worflow for variant calling with freebayes

// We need in term to be able to specify options ... as this 
/*
options.args = [ 
    "--db ${params.abricate_db}",
    "--minid ${params.minid}",
    "--mincov ${params.mincov}"
    ]
*/


include { BWAMEM_VAR } from "../modules/BWAMEM.nf"
include { SAMTOOLS_INDEX } from "../modules/SAMTOOLS.nf"


workflow MAPPING_VAR {

    take: 
    illumina_input_ch // channel: [ID, R1, R2]
    ref_ch // channel: assembly



    main:

    // TODO - test if better if we clean/trim reads or not So we need to have a read quality module

    BWAMEM_VAR(illumina_input_ch.combine(ref_ch)) 
    SAMTOOLS_INDEX(BWAMEM_VAR.out.shortbam_ch) 



    emit:
    bam           = SAMTOOLS_INDEX.out.shortindexbam_ch // tuple: [ID, path(_short_bam)]
    bam_bai       = SAMTOOLS_INDEX.out.bam_bai  // tuple: [ID, path(bam_bai)]



}

