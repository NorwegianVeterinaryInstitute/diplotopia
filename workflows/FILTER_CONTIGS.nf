// Local: Sub-workflows
include { FILTER_CONTIGS_SUB } from '../subworkflows/FILTER_CONTIGS_SUB.nf'

workflow FILTER_CONTIGS {

    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:['ID', 'assembly'], skip: 1, sep:",", strip:true)
        .map { row -> (ID, assembly) =  [row.ID, row.assembly]}

    FILTER_CONTIGS_SUB(input_ch)

    
}

