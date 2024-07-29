// (population) Variants At Regions Where Reference Reads Uniquely Map module
// subworkflows
include {VAR_PREP_REF} from '../subworkflows/VAR_PREP_REF.nf'
include {VAR_PREP_SAMPLES} from '../subworkflows/VAR_PREP_SAMPLES.nf'
include { CALL_SAMPLES; CALL_POP; COMBINE_INDIV_CALLS }  from "../subworkflows/VAR_FREEBAYES_SAMPLES_POP.nf"

// modules
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_SAMPLES;
            BCFTOOLS_FILTER as BCFTOOLS_FILTER_POP;
            BCFTOOLS_MERGE } from "../modules/BCFTOOLS.nf"

workflow VARWRRUM {

    if (params.phase == "prep_ref" | params.phase == "prep_samples" | params.phase == "call_variants" ) {

        // Part 0: Preparation: Reference - only illumina reads - required to set up parameters for mapping
        ref_input_ch = Channel.fromPath(params.ref_input, checkIfExists: true)
            .splitCsv(header:['ID', 'assembly', 'R1', 'R2', 'comments'], skip: 1, sep:",", strip:true)
            .map { row -> (ID, assembly, R1, R2) =  [row.ID, row.assembly, row.R1, row.R2 ]}

        VAR_PREP_REF(ref_input_ch)

    }


    // Part 0: Preparation : samples - only illunmina reads 
    // TODO make a if statement to check if previous done 


    if (params.phase == "prep_samples" | params.phase == "call_variants" ) {
        // TODO check if output from prep ref exist 
        // https://nextflow-io.github.io/patterns/process-when-empty/
        
        ifEmpty(VAR_PREP_REF.out.masked_ref) {
            error "Reference preparation not done yet"
        }

        input_ch = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header:['ID', 'R1', 'R2', 'ploidy', 'comments'], skip: 1, sep:",", strip:true)
            .map { row -> (ID, R1, R2, ploidy ) =  [row.ID, row.R1, row.R2 , row.ploidy]}

        illumina_input_ch = 
            input_ch
            .map { (ID, R1, R2) =  [it[0], it[1], it[2]] } 

        VAR_PREP_SAMPLES(illumina_input_ch, VAR_PREP_REF.out.masked_ref)
    }

    if (params.phase == "call_variants" ) {
        ifEmpty(VAR_PREP_SAMPLES.out.bam_bai) {
            error "Preparation samples not done yet"
        }

        // Calling variants 
        // --------------------------------
        // for samples individually 
        ploidy_ch = input_ch
                    .map{ (ID, ploidy) = [it[0], it[3]] }

        CALL_SAMPLES(
            VAR_PREP_SAMPLES.out.bam_bai, 
            ploidy_ch,  
            VAR_PREP_SAMPLES.out.ref_index
            )

        // --------------------------------
        // for population samples (from merged bam)
        // TODO now only ploidy 2
 
        CALL_POP(
            VAR_PREP_SAMPLES.out.pop_id, 
            VAR_PREP_SAMPLES.out.pop_bam, 
            VAR_PREP_SAMPLES.out.pop_bai, 
            VAR_PREP_SAMPLES.out.ref_index
            ) 

        // --------------------------------
        // combined individual vcf files (merged as population ) 
        COMBINE_INDIV_CALLS(
            CALL_SAMPLES.out.norm_vcfgz_index, 
            VAR_PREP_SAMPLES.out.ref_index.map{ (ref) = [it[0]] }
            )
                

    }
    
}





