include { MINIMAP } from "../modules/MINIMAP.nf"
include { TOBAM } from "../modules/TOBAM.nf"
include { QUALIMAP } from "../modules/QUALIMAP.nf"
include { MULTIQC } from "../modules/MULTIQC.nf"
include { QUAST_ASSEMBLY } from "../modules/QUAST.nf"
include { QUAST_REF } from "../modules/QUAST.nf"

workflow COVERAGE {
    
    //input csv format ID, path assembly, path reads
    if (!params.input) { exit 1, "Missing input file"}

    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:['ID', 'assembly', 'reads', 'comments'], skip: 1, sep:",", strip:true)
        .map { row -> (ID, assembly, reads) =  [ row.ID, row.assembly, row.reads ]}

    //input_ch.view()
    // for long reads only 
    // how to do by combining both short an long reads ? 
    //bwa mem
    // merge sam/bam ? 

    MINIMAP(input_ch)
    TOBAM(MINIMAP.out.sam_ch)
    QUALIMAP(TOBAM.out.bam_ch)
    MULTIQC(QUALIMAP.out.qualimap_ch)

    // need combining channel by ID to be able to keep cohesion
    quast_ch = input_ch.combine(TOBAM.out.bam_ch, by: 0)

    // if provide reference use ref otherwise use assembly as ref
    if(params.quast_ref){
        // not tested
        quast_ref_ch = Channel
            .fromPath(params.quast_ref, checkIfExists: true)
        quast_annot_ch = Channel
            .fromPath(params.quast_annot, checkIfExists: true)
        
        //quast_ref_ch.view()

        QUAST_REF(quast_ch, quast_ref_ch, quast_annot_ch)
            
            } // if no reference use assembly
            else {
                QUAST_ASSEMBLY(quast_ch)
            }

}
