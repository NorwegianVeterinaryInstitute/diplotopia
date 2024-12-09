/ subworkflows
include {} from '../subworkflows/VAR_PREP_REF.nf'
include {} from '../subworkflows/VAR_PREP_SAMPLES.nf'
include {}  from "../subworkflows/VAR_FREEBAYES_SAMPLES_POP.nf"

// modules
include { } from "../modules/BCFTOOLS.nf"

workflow FunGen {

    // input csv format ID, path assembly, path reads
    if (!params.input) { exit 1, "Missing input file"}

    // input required - short reads (illumina)
    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:['ID', 'R1', 'R2', 'comments'], skip: 1, sep:",", strip:true)
        .map { row -> (ID, R1, R2) =  [row.ID, row.R1, row.R2 ]}
 
        // ref genome given as param 

    if (params.phase == "prep_ref" | params.phase == "prep_samples" | params.phase == "call_variants" ) {

        // Part 0: Preparation: Reference - only illumina reads - required to set up parameters for mapping
        ref_input_ch = Channel.fromPath(params.ref_input, checkIfExists: true)
            .splitCsv(header:['ID', 'assembly', 'R1', 'R2', 'comments'], skip: 1, sep:",", strip:true)
            .map { row -> (ID, assembly, R1, R2) =  [row.ID, row.assembly, row.R1, row.R2 ]}

        VAR_PREP_REF(ref_input_ch)

    
    }

}

